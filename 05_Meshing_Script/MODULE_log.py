# ============================================================
# CONSOLE LOG STYLING
# ============================================================
# One consistent style for status output across the meshing pipeline,
# instead of ad-hoc print() formatting at each call site:
#   info()    - blue,  general status/progress
#   success() - green, a step or computed value that completed as expected
#   warn()    - red,   something the pipeline worked around and you should
#                       be aware of (not fatal — a fatal problem is a raised
#                       exception, not a warning)
# All three labels pad to the same width so the messages that follow them
# line up in a ragged-right column regardless of which one fires.

_RESET = "\033[0m"
_BOLD_BLUE = "\033[1;34m"
_BOLD_BRIGHT_BLUE = "\033[1;94m"
_BOLD_GREEN = "\033[1;32m"
_BOLD_RED = "\033[1;31m"


def info(message):
    print(f"{_BOLD_BLUE}INFO:{_RESET}    {message}")


def success(message):
    print(f"{_BOLD_GREEN}SUCCESS:{_RESET} {message}")


def warn(message):
    print(f"{_BOLD_RED}WARNING:{_RESET} {message}")


def subtitle(message):
    print(f"\n{_BOLD_BRIGHT_BLUE}{message}{_RESET}")
