#!/usr/bin/env python3
"""
pkl.py

Tools for constructing P(K)L sequences and verifying that sequences satisfy
P(K)L constraints.

MIT License

Copyright (c) 2021 by Abhinav Nellore and Rachel Ward

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from math import ceil, floor, log

def check(seq, low_mem=False, temp_dir=None):
    """ Checks if seq is a P(K)L sequence.

        seq: string or list of integers to check
        low_mem: True iff check should use sqlite dictionary on disk rather
            than memory
        temp: path to temporary directory or None if default should be used

        Return value: True iff seq is a P(K)L sequence.
    """
    L = len(seq)
    K = len(set(seq))
    if low_mem:
        from sqlitedict import SqliteDict
        from shutil import rmtree
        from tempfile import mkdtemp
        from os.path import join
        temp_dir = mkdtemp(dir=temp_dir)
    for n in range(1, ceil(log(L, K)) + 1):
        if low_mem:
            counts = SqliteDict(join(temp_dir, "pkl.sqlite"), autocommit=True)
        else:
            counts = {}
        for i in range(L):
            current_substring = tuple(seq[i:i + n] + seq[
                    :i + n - L if i + n - L > 0 else 0
                ])
            if current_substring not in counts:
                counts[current_substring] = 1
            else:
                counts[current_substring] += 1
        possibles = (ceil(L / K**n),
                     floor(L / K**n))
        if n == ceil(log(L, K)):
            if len(counts) != L:
                if low_mem:
                    rmtree(temp_dir)
                return False
        else:
            if len(counts) != K**n:
                if low_mem:
                    rmtree(temp_dir)
                return False
            for substring in counts:
                if counts[substring] not in possibles:
                    if low_mem:
                        rmtree(temp_dir)
                    return False
    if low_mem:
        rmtree(temp_dir)
    return True

if __name__ == "__main__":
    import argparse
    from sys import stdin
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
            help=("construct a P(K)L sequence or check that a sequence is a "
                  "P(K)L sequence"),
            dest='command'
        )
    check_parser = subparsers.add_parser("check")
    construct_parser = subparsers.add_parser("construct")
    check_parser.add_argument('-s', '--separator', type=str, metavar='<str>',
            default='',
            help=("separator between successive characters; to specify that "
                  "this is a newline in, e.g., Bash, use $'\n'")
        )
    check_parser.add_argument('-i', '--input', type=str, metavar='<file>',
            default=None,
            help=("file from which to read sequence; use '-' or leave"
                  "unspecified if stdin")
        )
    check_parser.add_argument('-l', '--low-memory', action='store_const',
            const=True,
            default=False,
            help=("use less memory when checking if sequence satisfies"
                  "P(K)L constraints")
        )
    check_parser.add_argument('-t', '--temp-dir', metavar="<dir>",
            default=None,
            help=("directory in which to store temporary files; defaults to "
                  "Python3's chain under tempfile.gettempdir() at "
                  "https://docs.python.org/3/library/tempfile.html")
        )
    construct_parser.add_argument('-k', '--alphabet-size', type=int,
            metavar='<int>',
            default=2,
            help="alphabet size")
    construct_parser.add_argument('-l', '--sequence-length', type=int,
            required=True,
            metavar='<int>',
            help='sequence length')
    args = parser.parse_args()
    if args.command == 'check':
        with open(args.input) if args.input != '-' else stdin as input_stream:
            seq = input_stream.read().strip().split(args.separator)
        if not check(seq, low_mem=args.low_memory, temp_dir=args.temp_dir):
            raise RuntimeError('Input is NOT a P(K)L sequence.')
        print("Input is a P(K)L sequence", file=stderr)

