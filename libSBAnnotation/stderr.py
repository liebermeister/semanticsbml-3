"""Default console feedback functions (messages, questions, progress)."""
import os
import sys

# set SEMANTICSBML_QUIET=1 to suppress progress messages
def _quiet():
    return os.environ.get('SEMANTICSBML_QUIET', '') not in ('', '0')

def info(message,close=False,clear=False):
    if _quiet(): return
    sys.stderr.write(message)
    sys.stderr.flush()

def yes_or_no(prompt = ''):
    # never block when not run from an interactive terminal
    if not sys.stdin or not sys.stdin.isatty():
        return False
    while 1:
        try:
            return input(prompt + ' (y/n) ').strip()[0] in ('y', 'Y')
        except IndexError:
            sys.stderr.write('Please answer y or n\n')
        except EOFError:
            return False

def downloadProgress(size_known,num_done,start=False,stop=False,url=''):
    """internal default function if no external UI functions were set"""
    if _quiet(): return
    if start:
        sys.stderr.write("Loading %s ..." % url)
    elif stop:
        sys.stderr.write(" done\n")
    sys.stderr.flush()

status_fkt = info
yes_or_no_fkt = yes_or_no
