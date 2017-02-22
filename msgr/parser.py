"""
parser.py
Parse the input code from execution program.
"""

import argparse
import sys


def parse(argv):
    '''
    Parse the input command line.

    Parameters
    ==========
        argv: str
            Input command line string.
    Notes
    =====

    Returns directory and field for initialization.
    '''

    parser = argparse.ArgumentParser(
        description="Start MSGR - volume Matching Satellite and Ground Radar.")

    parser.add_argument(
        '-s', '--script', type=str,
        help='Path to configuration file.', default=None)

    # Parse the args
    args = parser.parse_args(argv)

    if args.script is None:
        sys.exit()

    return args.script
