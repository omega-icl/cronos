#!/usr/bin/env python3
"""check_version.py -- CRONOS's version is written in four places, and the MC++ release it is built against in
four; each set must agree.

CRONOS's own version:
    pyproject.toml                  version = "X.Y.Z"           what pip and the wheels report
    src/interface/cronos.cpp        kVersion = "X.Y.Z"          what `cronos.__version__` reports
    doc/CRONOS.dox                  PROJECT_NUMBER              what the Doxygen pages are titled with
    src/makefile                    export version = X.Y.Z      what the makefile build prints

The MC++ release the module is built against (shared types cross between the two modules, so the headers the
wheels compile against and the pymcpp they run against must be the SAME release):
    .github/workflows/wheels.yml    MCPP_REF: vA.B.C            the MC++ sources the wheels compile against
    pyproject.toml                  pymcpp~=A.B.C               the runtime dependency
    pyproject.toml                  pymcpp==A.B.C  (twice)      the wheel builds' build requirements
    --mcpp-root <dir>               that tree's pyproject.toml  the MC++ a local build compiles against (optional)

The first CRONOS v5.1.0 tag had MCPP_REF at v5.1.0 and the three pins at 5.0.4, and `check` passed: the wheel
compiled against 5.1.0 headers, installed pymcpp 5.0.4, and failed its own import check in CI (2026-10-27).
Run from the repository root, or with the root as the first argument; exits non-zero on any disagreement,
naming each file's value.  With --expect X.Y.Z, also requires that CRONOS version (the release workflow passes
the tag's).
"""
import os
import re
import sys


def grab(path, pattern, flags=re.M):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        m = re.search(pattern, f.read(), flags)
    return m.group(1) if m else None


def grab_all(path, pattern, flags=re.M):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as f:
        return re.findall(pattern, f.read(), flags)


def main(argv):
    root, expect, mcpp_root = ".", None, None
    args = list(argv[1:])
    while args:
        a = args.pop(0)
        if a == "--expect":
            expect = args.pop(0)
        elif a == "--mcpp-root":
            mcpp_root = args.pop(0)
        else:
            root = a
    rc = 0

    # ---- CRONOS's own version
    own = [
        ("pyproject.toml version", grab(os.path.join(root, "pyproject.toml"), r'^version\s*=\s*"([^"]+)"')),
        ("cronos.cpp kVersion", grab(os.path.join(root, "src", "interface", "cronos.cpp"), r'kVersion\s*=\s*"([^"]+)"')),
        ("doc/CRONOS.dox PROJECT_NUMBER", grab(os.path.join(root, "doc", "CRONOS.dox"), r'^PROJECT_NUMBER\s*=\s*(\S+)')),
        ("src/makefile version", grab(os.path.join(root, "src", "makefile"), r'^export\s+version\s*=\s*(\S+)')),
    ]
    print("CRONOS version:")
    for name, v in own:
        print("  %-34s %s" % (name, v if v else "(not found)"))
    values = {v for _, v in own}
    if len(values) != 1 or None in values:
        print("check_version.py: FAIL -- CRONOS's version is not the same everywhere"); rc = 1
    elif expect is not None and expect not in values:
        print("check_version.py: FAIL -- the files say %s, the tag says %s" % (values.pop(), expect)); rc = 1
    else:
        print("  the CRONOS version is %s everywhere" % next(iter(values)))

    # ---- the MC++ release
    pins = [
        ("wheels.yml MCPP_REF", grab(os.path.join(root, ".github", "workflows", "wheels.yml"), r'^\s*MCPP_REF:\s*v?([0-9]+\.[0-9]+\.[0-9]+)')),
        ("pyproject.toml pymcpp~=", grab(os.path.join(root, "pyproject.toml"), r'"pymcpp~=([0-9]+\.[0-9]+\.[0-9]+)"')),
    ]
    for k, v in enumerate(grab_all(os.path.join(root, "pyproject.toml"), r'"pymcpp==([0-9]+\.[0-9]+\.[0-9]+)"')):
        pins.append(("pyproject.toml pymcpp== (%d)" % (k + 1), v))
    if mcpp_root:
        pins.append(("MCPP_ROOT pyproject.toml version", grab(os.path.join(mcpp_root, "pyproject.toml"), r'^version\s*=\s*"([^"]+)"')))
    print("MC++ release built against:")
    for name, v in pins:
        print("  %-34s %s" % (name, v if v else "(not found)"))
    pvals = {v for _, v in pins}
    if len(pvals) != 1 or None in pvals:
        print("check_version.py: FAIL -- the MC++ release is not the same everywhere: the wheels compile against MCPP_REF"
              " and run against the pymcpp pin, and they must be one release"); rc = 1
    else:
        print("  the MC++ release is %s everywhere" % next(iter(pvals)))
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))
