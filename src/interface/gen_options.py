#!/usr/bin/env python3
# gen_options.py -- generate pybind11 bindings for a CRONOS Options struct from its header.
#
#   gen_options.py HEADER STRUCT_LINE|auto CXX_TYPE FUNC           one struct, to stdout
#   gen_options.py --all  [SRC_DIR [OUT_DIR]]                        the four option structs of CRONOS, written
#   gen_options.py --check [SRC_DIR [OUT_DIR]]                       ... compared with the committed files (exit 1 if stale)
# STRUCT_LINE `auto` finds the struct by name (the first `struct Options`); the generated header then carries no line
# number, so unrelated edits of a header do not make its generated file stale.  SRC_DIR defaults to the script's parent
# (src/), OUT_DIR to the script's own directory (src/interface/).  --all/--check format with clang-format when it is
# found (MC++'s .clang-format, via MCPP_ROOT or the CRONOS_CLANG_FORMAT_STYLE variable); --check is skipped without it.
#
# Parses the struct starting at STRUCT_LINE (1-based, the "struct Options" line) of HEADER:
#   * data members  "TYPE NAME;"  (initialisers ignored), docs from trailing //!< (and //!< continuation lines) or
#     the preceding //! block;
#   * nested structs "struct t_X { ... } NAME;"  -> a nested Python class t_X and a read-write attribute NAME;
#   * enums declared in the struct "enum E { A = 0, //!< doc ... };"  -> py::native_enum on the Options class;
#   * aliases "using E = Q;" + "static constexpr E V = Q::W;"  -> a native_enum of Q with the re-exported values.
# Static members, methods and constructors are skipped.  Emits one C++ function
#   template <class PyOpt> void FUNC( PyOpt& c )
# to be called on the py::class_ of CXX_TYPE.  Field types that are enums defined elsewhere (e.g. in the base's
# Options) are bound where they are declared; the report on stderr lists every field and flags unknown types.
import os, re, sys

def strip_comment(l): return re.sub(r'//.*', '', l)

def parse_block(lines, i0):
    """lines[i0] contains the opening '{' of a struct; return (items, index after the closing '};' / '} NAME;')."""
    depth = 0; i = i0; items = []; doc = []
    while True:
        l = lines[i]; s = strip_comment(l)
        if depth == 1 and l.strip().startswith('#'):
            items.append(('pp', l.strip(), None, [], None)); i += 1; continue
        if depth == 1:
            st = l.strip()
            if st.startswith('//!') and not st.startswith('//!<'):
                t = re.sub(r'^//!\s?(@brief\s*)?', '', st); doc.append(t)
            elif st.startswith('//!<') and items and items[-1][0] == 'field':
                items[-1][3].append(re.sub(r'^//!<\s?', '', st))
            m = re.match(r'^\s*struct\s+(\w+)\s*$', l) or re.match(r'^\s*struct\s+(\w+)\s*\{', l)
            if m:
                j = i if '{' in l else i + 1
                sub, k = parse_block(lines, j)
                mm = re.match(r'^\s*\}\s*(\w+)\s*;', lines[k - 1])
                items.append(('struct', m.group(1), mm.group(1) if mm else None, doc, sub)); doc = []
                i = k; continue
            m = re.match(r'^\s*enum\s+(class\s+)?(\w+)', l)
            if m:
                vals = []; j = i
                while '}' not in strip_comment(lines[j]) or j == i and '{' not in lines[j]:
                    j += 1
                    if lines[j].strip().startswith('#'):
                        vals.append(('#', lines[j].strip())); continue
                    for mv in re.finditer(r'\b([A-Z][A-Z_0-9]*)\s*(=\s*[-\w]+)?\s*,?\s*(//!<\s*(.*))?$', lines[j].strip() if j <= len(lines) else ''):
                        if strip_comment(lines[j]).strip() and not strip_comment(lines[j]).strip().startswith('}'):
                            vals.append((mv.group(1), (mv.group(4) or '').strip()))
                items.append(('enum', m.group(2), bool(m.group(1)), doc, vals)); doc = []
                i = j + 1; continue
            m = re.match(r'^\s*using\s+(\w+)\s*=\s*([\w:]+)\s*;', s)
            if m: items.append(('alias', m.group(1), m.group(2), doc, [])); doc = []; i += 1; continue
            m = re.match(r'^\s*static\s+constexpr\s+(\w+)\s+(\w+)\s*=\s*([\w:]+)\s*;', s)
            if m:
                for it in items:
                    if it[0] == 'alias' and it[1] == m.group(1): it[4].append((m.group(2), m.group(3)))
                i += 1; continue
            m = re.match(r'^\s*(?!return|static|typedef|using|friend|virtual|explicit|inline)([A-Za-z_][\w:<>, ]*?[\w>])\s+([A-Za-z_]\w*)\s*(=[^;]*)?;\s*(//!<\s*(.*))?$', l)
            if m and '(' not in s:
                items.append(('field', m.group(2), m.group(1).strip(), doc + ([m.group(5).strip()] if m.group(5) else []))); doc = []
                i += 1; continue
            if st and not st.startswith('//'): doc = []
        depth += s.count('{') - s.count('}')
        i += 1
        if depth == 0 and '}' in s: return items, i

def doc_str(d):
    t = ' '.join(x for x in d if x).strip()
    t = re.sub(r'\s+', ' ', t).replace('\\', '\\\\').replace('"', '\\"')
    return t or 'undocumented'

def emit(items, cxx, pyvar, out, report, prefix, enums_known):
    for it in items:
        if it[0] == 'enum':
            _, name, is_class, doc, vals = it
            out.append('  py::enum_<%s::%s>( %s, "%s", "%s" )' % (cxx, name, pyvar, name, doc_str(doc)))
            for v, vd in vals:
                if v == '#': out.append(vd); continue
                out.append('    .value( "%s", %s::%s%s, "%s" )' % (v, cxx, (name + '::') if is_class else '', v, doc_str([vd]) if vd else v))
            out.append('    %s;' % ('' if is_class else '.export_values()'))
            enums_known.add(name)
        elif it[0] == 'alias' and it[4]:
            _, name, target, doc, vals = it
            target = target if target.startswith('mc::') else 'mc::' + target
            out.append('  py::enum_<%s>( %s, "%s", "%s" )' % (target, pyvar, name, doc_str(doc) if doc != [] else name + ' (' + target + ')'))
            for v, q in vals: out.append('    .value( "%s", %s, "%s" )' % (v, q if q.startswith('mc::') else 'mc::' + q, q))
            out.append('    .export_values();')
            enums_known.add(name)
    for it in items:
        if it[0] == 'pp':
            out.append(it[1]); continue
        if it[0] == 'field':
            _, name, typ, doc = it
            out.append('  %s.def_readwrite( "%s", &%s::%s, "%s" );' % (pyvar, name, cxx, name, doc_str(doc)))
            report.append('%-40s %-26s %s' % (prefix + name, typ, '' if typ in ('bool','int','double','unsigned','unsigned int','size_t','std::string','long','sunrealtype','float') or typ in enums_known or typ.split('::')[-1] in enums_known else '  <-- type bound elsewhere?'))
        elif it[0] == 'struct' and it[2]:
            _, tname, fname, doc, sub = it
            v = 'c_' + fname
            out.append('  py::class_<%s::%s> %s( %s, "%s", "%s" );' % (cxx, tname, v, pyvar, tname, doc_str(doc)))
            out.append('  %s.def( py::init<>() );' % v)
            emit(sub, '%s::%s' % (cxx, tname), v, out, report, prefix + fname + '.', enums_known)
            out.append('  %s.def_readwrite( "%s", &%s::%s, "%s" );' % (pyvar, fname, cxx, fname, doc_str(doc)))

def generate( hdr, line, cxx, func ):
    lines = open(hdr).read().split('\n')
    if line == 'auto':
        i = next( k for k, l in enumerate(lines) if re.match(r'^\s*struct Options\b', l) )
    else:
        i = int(line) - 1
    while '{' not in lines[i]: i += 1
    items, _ = parse_block(lines, i)
    out = ['// Copyright (C) Benoit Chachuat, Imperial College London.', '// All Rights Reserved.',
           '// This code is published under the Eclipse Public License.', '',
           '// GENERATED by gen_options.py from %s -- do not edit; regenerate' % hdr.split('/')[-1],
           '// (cmake --build <dir> --target gen-options).' if line == 'auto' else '// (struct at line %d).' % int(line),
           '', '#pragma once', '', '#include <pybind11/pybind11.h>', '', 'namespace py = pybind11;', '',
           'template <class PyOpt>', 'void', '%s(PyOpt& c)' % func, '{']
    report = []; known = set()
    emit(items, cxx, 'c', out, report, '', known)
    out.append('}')
    return '\n'.join(out) + '\n', report

SPECS = [ ('ffmodel.hpp',       'mc::FFModel::Options',        'bind_ffmodel_options', 'gen_ffmodel_options.hpp'),
          ('ocfeslv.hpp',       'mc::OCFESLV::Options',        'bind_ocfeslv_options', 'gen_ocfeslv_options.hpp'),
          ('odeslv_cvodes.hpp', 'mc::ODESLV_CVODES::Options',  'bind_odeslv_options',  'gen_odeslv_options.hpp'),
          ('base_cvodes.hpp',   'mc::BASE_CVODES::Options',    'bind_cvodes_options',  'gen_cvodes_options.hpp') ]

def clang_format( text ):
    import shutil, subprocess
    exe = shutil.which('clang-format')
    style = os.environ.get('CRONOS_CLANG_FORMAT_STYLE')
    if not style and os.environ.get('MCPP_ROOT') and os.path.exists(os.environ['MCPP_ROOT'] + '/.clang-format'):
        style = 'file:' + os.environ['MCPP_ROOT'] + '/.clang-format'
    if not exe or not style: return None
    r = subprocess.run([exe, '-style=' + style, '--assume-filename=x.hpp'], input=text, capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None

def main():
    a = sys.argv[1:]
    if a and a[0] in ('--all', '--check'):
        here = os.path.dirname(os.path.abspath(__file__))
        src = a[1] if len(a) > 1 else os.path.join(here, '..')
        dst = a[2] if len(a) > 2 else here
        stale = []
        for hpp, cxx, func, outname in SPECS:
            text, _ = generate(os.path.join(src, hpp), 'auto', cxx, func)
            fmt = clang_format(text)
            if fmt is None:
                if a[0] == '--check':
                    print('gen_options.py --check: skipped (clang-format or its style not found: set MCPP_ROOT)'); return 0
                sys.stderr.write('gen_options.py: clang-format not found -- %s written UNFORMATTED\n' % outname); fmt = text
            path = os.path.join(dst, outname)
            if a[0] == '--all':
                open(path, 'w').write(fmt); print('wrote ' + path)
            elif not os.path.exists(path) or open(path).read() != fmt:
                stale.append(outname)
        if a[0] == '--check':
            if stale:
                print('gen_options.py --check: STALE (regenerate with the gen-options target): ' + ', '.join(stale)); return 1
            print('gen_options.py --check: the %d generated option headers are up to date' % len(SPECS))
        return 0
    text, report = generate(a[0], a[1], a[2], a[3])
    print(text, end='')
    sys.stderr.write('%s: %d fields\n' % (a[2], len(report)) + '\n'.join('  ' + r for r in report) + '\n')
    return 0

sys.exit(main())
