#!/usr/bin/env python3
# Copyright (C) Benoit Chachuat, Imperial College London.
# All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# gen_environment.py -- docs/ENVIRONMENT.md, GENERATED from the code and a curated registry (2026-10-07, WORKPLAN 3.D).
#
#   gen_environment.py --write [SRC_DIR [DOC]]    write the document
#   gen_environment.py --check [SRC_DIR [DOC]]    exit 1 if a header reads a CRONOS_* variable the registry does not
#                                                  document, if the registry documents one no header reads, or if the
#                                                  document is not what --write would produce
#
# SRC_DIR defaults to the script's parent (src/), DOC to ../../docs/ENVIRONMENT.md relative to the script.  A new
# environment variable therefore cannot be added silently: `make check` fails until it is described here.  Read sites are
# recorded by HEADER, not by line, so that ordinary edits of a header do not make the document stale.
import os, re, sys

# ---- the registry: variable -> (group, what it does).  Written from the code (the comment at each read site). ----------
GROUPS = [
  ( 'resource', 'Resources',
    'Sets the default of an option that does not change results.' ),
  ( 'audit', 'Setup reports and audits',
    'Print reports at setup; none changes what is computed.  Those marked "default of X" set the default of an option, '
    'which a driver setting the option overrides.' ),
  ( 'plan', 'Interface-plan dumps', 'Print the interface plan\'s decisions, at display level 3.' ),
  ( 'rank', 'Rank and factorisation diagnostics', 'Cross-check or report the numerical rank decisions; print only.' ),
  ( 'solve', 'Solve diagnostics', 'Report on solutions and sensitivities; print only.' ),
  ( 'test', 'Test switches (not API)',
    'Change behaviour; kept ONLY for the corpus drivers that exercise them.  A model should not set them.' ),
]
REGISTRY = {
  # resources
  'CRONOS_MAXTHREAD': ( 'resource', 'Default of `MAXTHREAD` (0 = no cap): lets a cluster job cap the threads without editing code.' ),
  # setup reports and audits
  'CRONOS_AUDIT': ( 'audit', 'Umbrella for the cheap setup reports: default of `AUDIT.EQN_LEGEND`, `AUDIT.COVERAGE` and `AUDIT.RANK_SYMMETRY` (not the spectrum, which runs a dense SVD).' ),
  'CRONOS_AUDIT_LEGEND': ( 'audit', 'Default of `AUDIT.EQN_LEGEND`: the equation legend.' ),
  'CRONOS_AUDIT_COVERAGE': ( 'audit', 'Default of `AUDIT.COVERAGE`: the coverage report.' ),
  'CRONOS_AUDIT_RANK': ( 'audit', 'Default of `AUDIT.RANK_SYMMETRY`: the rank-symmetry report.' ),
  'CRONOS_AUDIT_SPECTRUM': ( 'audit', 'Default of `AUDIT.SPECTRUM`: the spectrum report (a dense SVD; not implied by `CRONOS_AUDIT`).' ),
  'CRONOS_SPECTRUM_MAX': ( 'audit', 'Size cap of the spectrum report (`AUDIT.SPECTRUM_MAX`).' ),
  'CRONOS_AUDIT_DETERMINACY': ( 'audit', 'Default of `DETERMINACY.AUDIT`: the numeric (SPQR) determinacy audit, in place of the structural determinacy line.' ),
  'CRONOS_AUDIT_MAX_N': ( 'audit', 'Budget of the determinacy audit: above this many rows or columns it reports SKIPPED instead of factorising (default unlimited).' ),
  'CRONOS_DETERMINACY_GAP_MIN': ( 'audit', 'Determinacy audit: a gap at the rank cut below this triggers a re-evaluation at 1.3x the reference point (default 1e2).' ),
  'CRONOS_KPRIME_STRUCTURAL': ( 'audit', 'With the numeric determinacy audit on, print the structural determinacy line beside it (a corpus-wide equality check).' ),
  'CRONOS_AUDIT_SETUP_ONLY': ( 'audit', '`setup()` returns false right after the audit reports, before the derivative caches and any solve -- for structural corpus sweeps.' ),
  'CRONOS_AUDIT_TAG': ( 'audit', 'A one-line machine-readable `[cov1]` summary per setup, for corpus sweeps.' ),
  'CRONOS_HEADER_FULL': ( 'audit', 'Print the solver header\'s revision summary at setup (display level 2).' ),
  'CRONOS_AUDIT_SYMROW': ( 'audit', 'Model layer: report the value-slaved states of the principal-symbol row pass.' ),
  'CRONOS_AUDIT_EVOFACE': ( 'audit', 'Evolution-face probes: the hoisted evolution-domain detection and the claims on evolution faces.' ),
  'CRONOS_AUDIT_IFACEPLAN': ( 'audit', 'Print the interface-plan report at setup, without a driver edit.' ),
  'CRONOS_AUDIT_DEPDUMP': ( 'audit', 'Dump each equation\'s dependency map (key, resolved name, type).' ),
  'CRONOS_AUDIT_RECVSCAN': ( 'audit', 'Static scan at setup: for each claim, whether a consuming row carries a face-direction derivative.' ),
  'CRONOS_AUDIT_AUXRECV': ( 'audit', 'Report the auxiliary-receiver test of the weak-SAT receiver states.' ),
  'CRONOS_AUDIT_SITESTATES': ( 'audit', 'Report the states at each claim site, with the reduced-flux continuity predicate.' ),
  'CRONOS_AUDIT_MINT': ( 'audit', 'Report why candidate continuity claims were not minted (e.g. no continuity row could be built).' ),
  'CRONOS_AUDIT_TAUCOL': ( 'audit', 'The `[taucol/plan]` report: the plan\'s tau columns.' ),
  'CRONOS_AUDIT_DROPSET': ( 'audit', 'Name the claims the redundancy pass drops, beside the rank deficiency they address.' ),
  'CRONOS_AUDIT_DROPNULL': ( 'audit', 'Report the W-test of dropped claims: how many must be restored (rank W) and which (the pivot rows).' ),
  'CRONOS_KEEP_EXPLICIT_AUDIT': ( 'audit', 'Check that every explicit tau column arrives with its row (the keep-explicit invariant).' ),
  'CRONOS_REF_AUDIT': ( 'audit', 'Report the reference point each audit assembles at (setup, or the W-decide re-derivation).' ),
  # interface-plan dumps
  'CRONOS_PLAN_DUMP_STATE': ( 'plan', '`=<state>`: dump the plan\'s decisions for one state.' ),
  'CRONOS_PLAN_DUMP_CLAIMS': ( 'plan', '`=<state>`: claim census for one state -- which seam nodes carry a claim, and its fate (dropped, unrealised, forest-redundant).' ),
  'CRONOS_PLAN_RECV_CENSUS': ( 'plan', 'Receiver census by row role, per (state, direction) of the kept claims; for an auxiliary, how many receivers are its own LINK row.' ),
  'CRONOS_PLAN_DUMP_DETCOLS': ( 'plan', 'Report each determinacy-column invocation of the plan builder (its columns).' ),
  'CRONOS_PLAN_DUMP_DROPPED': ( 'plan', 'List the claims the plan drops (block, state, direction, element, face, kind).' ),
  # rank and factorisation
  'CRONOS_SPQR_STAT': ( 'rank', 'SPQR\'s statistics per factorisation: nnz(R), rank estimate, singletons, ordering used.' ),
  'CRONOS_DEBUG_RANK': ( 'rank', 'Umbrella: `CRONOS_SPQR_STAT` and `CRONOS_KEEP_EXPLICIT_AUDIT` together.' ),
  'CRONOS_STAIRCASE_VERIFY': ( 'rank', 'Run both live-column scans of the rank refinement (the fast one and its fallback) and report any disagreement.' ),
  'CRONOS_RANK_ORACLE': ( 'rank', 'Cross-check each numeric rank against a dense SVD of the same matrix at the same tolerance: both ranks and the singular values at the cut.' ),
  'CRONOS_RANK_ORACLE_MAXN': ( 'rank', 'Size cap of the rank oracle (default 2000).' ),
  'CRONOS_TAUNULL_CENSUS': ( 'rank', 'Name the collective null directions of the tau block when k\' > 0.' ),
  'CRONOS_STATENULL_CENSUS': ( 'rank', 'Name the state-side null directions (def_right > k\').' ),
  'CRONOS_STATENULL_ROWS': ( 'rank', 'With the state-null census: the rows that see each null direction.' ),
  'CRONOS_DUMP_COLS': ( 'rank', '`=c1,c2,...`: dump the nonzeros of selected columns of the audited Jacobian.' ),
  'CRONOS_DUMP_J': ( 'rank', '`=<file>`: write the audited Jacobian in coordinate form.' ),
  # solve diagnostics
  'CRONOS_SOLVE_VERBOSE': ( 'solve', 'Force `SOLVE.VERBOSE` on in every driver (force-on only), to trace a whole sweep\'s solves without editing drivers.' ),
  'CRONOS_AUDIT_BETA': ( 'solve', 'After each IC_TRACE solve with the projection applied: beta = Z^T lambda and the violation it absorbs.' ),
  'CRONOS_AUDIT_SENS': ( 'solve', '||dF/dvar|| restricted to the tau columns, once per sensitivity setup.' ),
  'CRONOS_DUP_SPREAD': ( 'solve', 'After every solve, per-state spreads over duplicate nodes: the interface conditions the solution actually satisfies.' ),
  'CRONOS_DUP_SPREAD_WARN': ( 'solve', 'Warning threshold of those spreads (default 1e-10; a report, not a failure).' ),
  # test switches
  'CRONOS_FORCE_BAD_KEEP_EXPLICIT': ( 'test', 'Fault injector for `OCFE_hybrid_adequacy`: makes the tau elimination stall, which is that driver\'s PASS condition.' ),
  'CRONOS_RESCUE_C2': ( 'test', 'Receiver policy of the C2/C2b study (default on; 0 forces the original +1 everywhere).  Kept for `OCFE_receiver` until the receiver-policy study concludes.' ),
  'CRONOS_WEAK_NATURAL_PENALTY': ( 'test', 'IC_WEAK receiver policy for continuity claims: 0, 1 or 2 (default).  Kept for `OCFE_receiver`, as above.' ),
}
MCPP = {   # read by MC++, not by CRONOS's headers -- documented, not scanned
  'CRONOS_FOLD_PARTIALS': 'MC++ (`ocbase.hpp`): default of the folding of chained partial derivatives (on); settable in code -- set it before building expressions.',
}

def scan( src ):
    found = {}
    for f in sorted( os.listdir( src ) ):
        if not f.endswith( '.hpp' ): continue
        for line in open( os.path.join( src, f ), errors='ignore' ):
            if line.strip().startswith( '//' ): continue
            for m in re.finditer( r'"(CRONOS_[A-Z0-9_]+)"', line ):
                found.setdefault( m.group( 1 ), set() ).add( f )
    return found

def render( found ):
    out = [ '# Environment variables', '',
            '<!-- GENERATED by src/interface/gen_environment.py -- do not edit; regenerate with the gen-environment target. -->', '',
            'CRONOS is configured through its **options** (`FFModel::Options`, `OCFESLV::Options`, the ODESLV options), which',
            'travel with the source.  The environment variables below are what remains (2026-10-07): diagnostics that instrument a',
            'whole corpus sweep without editing any driver, one resource setting, and three test switches kept for the drivers that',
            'exercise them.  **None of them is needed to run a model**, and none changes a result except the test switches.', '',
            'Most switches read `0`, `false`, `off` and `no` as off; a few are on as soon as they are SET, whatever the value --',
            'leave a variable unset to keep it off.  Variables shown with `=<...>` take a value.', '' ]
    for key, title, intro in GROUPS:
        rows = sorted( v for v, ( g, _ ) in REGISTRY.items() if g == key )
        out += [ '## ' + title, '', intro, '', '| variable | what it does | read in |', '|---|---|---|' ]
        for v in rows:
            out.append( '| `%s` | %s | %s |' % ( v, REGISTRY[ v ][ 1 ], ', '.join( '`%s`' % f for f in sorted( found.get( v, [] ) ) ) ) )
        out.append( '' )
    out += [ '## Read by MC++', '', '| variable | what it does |', '|---|---|' ]
    for v in sorted( MCPP ): out.append( '| `%s` | %s |' % ( v, MCPP[ v ] ) )
    out.append( '' )
    return '\n'.join( out )

def main():
    a = sys.argv[ 1: ]
    if not a or a[ 0 ] not in ( '--write', '--check' ):
        sys.stderr.write( 'usage: gen_environment.py --write|--check [SRC_DIR [DOC]]\n' ); return 2
    here = os.path.dirname( os.path.abspath( __file__ ) )
    src = a[ 1 ] if len( a ) > 1 else os.path.join( here, '..' )
    doc = a[ 2 ] if len( a ) > 2 else os.path.join( here, '..', '..', 'docs', 'ENVIRONMENT.md' )
    found = scan( src )
    undocumented = sorted( set( found ) - set( REGISTRY ) - set( MCPP ) )   # an MC++ header in SRC_DIR is documented below
    stale = sorted( set( REGISTRY ) - set( found ) )
    text = render( found )
    if a[ 0 ] == '--write':
        if undocumented or stale:
            print( 'gen_environment.py: undocumented %s, stale %s -- fix the registry first' % ( undocumented, stale ) ); return 1
        open( doc, 'w' ).write( text ); print( 'wrote %s (%d variables)' % ( doc, len( REGISTRY ) + len( MCPP ) ) ); return 0
    bad = []
    if undocumented: bad.append( 'read by a header but NOT documented: ' + ', '.join( undocumented ) )
    if stale: bad.append( 'documented but read by no header: ' + ', '.join( stale ) )
    if not bad and ( not os.path.exists( doc ) or open( doc ).read() != text ):
        bad.append( 'docs/ENVIRONMENT.md is out of date (regenerate with the gen-environment target)' )
    if bad:
        print( 'gen_environment.py --check: ' + '; '.join( bad ) ); return 1
    print( 'gen_environment.py --check: the %d environment variables are documented and docs/ENVIRONMENT.md is up to date' % len( REGISTRY ) )
    return 0

sys.exit( main() )
