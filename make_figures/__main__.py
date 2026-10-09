import sys

from . import main

argv = sys.argv[1:]
only = None
if "--only" in argv:
    i = argv.index("--only")
    only = {int(a) for a in argv[i + 1:] if a.isdigit()}
main(only=only, cpu="--cpu" in argv, fast="--fast" in argv)
