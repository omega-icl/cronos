#!/usr/bin/env python3
# Turn pyproject.toml (the EPL build of cronos-mcpp) into the GPL build's, in place -- CI only; never commit the result.
#   python tools/wheel_variant.py gpl     -> licence GPL-2.0-or-later, the GPL text shipped, UMFPACK and SPQR on
# The project NAME is unchanged: the GPL wheels are published as GitHub Release assets with a build tag (`wheel tags
# --build 1`), never on PyPI, and pip/uv prefer them when the release page is given (`pip install cronos-mcpp -f
# <release page>`).  Every edit must apply exactly once and the result is re-read with tomllib, so a pyproject.toml this
# script no longer matches -- or one already converted -- fails the build instead of producing a mislabelled wheel.
import re, sys, tomllib
def main():
    if sys.argv[1:] != [ "gpl" ]:
        sys.exit( "usage: wheel_variant.py gpl" )
    p = "pyproject.toml"; s = open( p ).read()
    name = tomllib.loads( s )[ "project" ][ "name" ]
    edits = [
        ( r'^license = "EPL-2\.0"',                       'license = "GPL-2.0-or-later"' ),
        ( r'^license-files = \[[^\]]*\]',                 'license-files = ["LICENSE", "LICENSES/GPL-2.0-or-later.txt", "THIRD_PARTY_NOTICES.txt"]' ),
        ( r'^description = "([^"]+)"',                    r'description = "\1 (GPL build, with UMFPACK and SPQR)"' ),
        ( r'^CRONOS_WITH_UMFPACK = "OFF"',                'CRONOS_WITH_UMFPACK = "ON"' ),
        ( r'^CRONOS_WITH_SPQR = "OFF"',                   'CRONOS_WITH_SPQR = "ON"' ),
    ]
    for pat, rep in edits:
        s, n = re.subn( pat, rep, s, flags = re.M )
        if n != 1: sys.exit( "wheel_variant.py: '%s' matched %d time(s), expected 1" % ( pat, n ) )
    t = tomllib.loads( s )
    d = t[ "tool" ][ "scikit-build" ][ "cmake" ][ "define" ]
    assert t[ "project" ][ "name" ] == name and t[ "project" ][ "license" ] == "GPL-2.0-or-later"
    assert d[ "CRONOS_WITH_UMFPACK" ] == "ON" and d[ "CRONOS_WITH_SPQR" ] == "ON"
    open( p, "w" ).write( s )
    print( "pyproject.toml -> the GPL build of %s (GPL-2.0-or-later; UMFPACK, SPQR on)" % name )
main()
