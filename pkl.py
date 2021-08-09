#!/usr/bin/env python3
"""
pkl.py

Tools for constructing P(K)L-sequences as described in

A. Nellore and R. Ward. "Arbitrary-length analogs to de Bruijn sequences"

and verifying that constructed sequences are P(K)L-sequences.

The implementation here is didactic. While it is uses O(L) space and takes
O(L log K) time, as discussed in the paper above, there are several obvious
optimizations that can be made, notably in searching for substrings at which to
join strings.

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
from os import remove
from fractions import gcd
from collections import deque

def check(seq, low_mem=False, temp_dir=None):
    """ Checks if a (circular) string is a P(K)L-sequence.

        seq: string or list of integers to check
        low_mem: True iff check should use sqlite dictionary on disk rather
            than memory
        temp: path to temporary directory or None if default should be used

        Return value: True iff seq is a P(K)L-sequence.
    """
    L = len(seq)
    K = len(set(seq))
    if low_mem:
        try:
            from sqlitedict import SqliteDict
        except ImportError:
            raise ImportError(
                "--low-memory/-l requires SqliteDict. Install it with "
                "'pip3 install -U sqlitedict'."
            )
        from shutil import rmtree
        from tempfile import mkdtemp
        from os.path import join
        from pickle import dumps
        from os import remove
        dumps_wrapper = lambda x: dumps(x, 0).decode("utf-8")
        temp_dir = mkdtemp(dir=temp_dir)
    else:
        dumps_wrapper = lambda x: x
    for n in range(1, ceil(log(L, K)) + 1):
        if low_mem:
            try:
                remove(join(temp_dir, "pkl.sqlite"))
            except FileNotFoundError:
                pass
            counts = SqliteDict(join(temp_dir, "pkl.sqlite"), autocommit=True)
        else:
            counts = {}
        for i in range(L):
            current_substring = dumps_wrapper(tuple(seq[i:i + n] + seq[
                    :i + n - L if i + n - L > 0 else 0
                ]))
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
            counts.close()
    if low_mem:
        rmtree(temp_dir)
    return True

def length_in_base(L, K):
    """ Converts integer to base-K number.

        L: length to convert
        K: desired base

        Return value: list of digits, from most to least significant
    """
    assert L >= 0
    digits = []
    while L:
        digits.append(L % K)
        L = L // K
    return digits[::-1]

def lempels_lift(seq, K):
    """ Constructs Lempel's lift of a sequence

        seq: list of integers to lift
        K: alphabet size

        Return value: list of lists of integers in Lempel's lift
    """
    mod_sum = 0
    L = len(seq)
    for i in range(L):
        mod_sum = (seq[i] + mod_sum) % K
    repetitions = K // gcd(mod_sum, K)
    lift = []
    for j in range(K * L // repetitions):
        lift.append(deque(0 for _ in range(repetitions * L)))
        lift[-1][-1] = j
        for i in range(L * repetitions):
            lift[-1][i] = (lift[-1][i-1] + seq[i % L]) % K
    return lift

def easy_join(lift, streak_length):
    """ Joins strings in lift on (x, x + 1, ..., x + streak_length - 1)

        Warning: string of form (x, x + 1, ..., x + streak_length - 1) and
        length streak_length should be unique to each string in lift

        streak_length: length of join string
        lift: strings to join, where the nth string is 1 + the (n-1)th string
    """
    string_length = len(lift[0])
    total_strings = len(lift)
    total_length = string_length * total_strings
    streak_length = ceil(log(total_length, K))
    streak = 0
    # So streaks can be detected starting at first character
    for i in range(string_length - streak_length + 2, string_length):
        if lift[0][
                (i - current) % string_length
            ] % K == (lift[0][
                    (i - current - 1) % string_length
                ] - 1) % K:
            streak += 1
        else:
            streak = 0
    current = 0
    seq = deque()
    # Join Lempel's lift shift-rule style
    for i in range(total_length):
        if lift[current % total_strings][
                    (i - current) % string_length
                ] % K == (
                lift[current % total_strings][
                    (i - current - 1) % string_length
                ] - 1) % K:
            streak += 1
        else:
            streak = 0
        seq.append(lift[current % total_strings][i - current])
        if streak == streak_length:
            # String of form (x, x+1, ..., x+n) found; continue on next
            # string
            current += 1
            # Streak should be reset
            streak = 0
    return seq

def hard_join()

def joined_lift(lift, digit):
    """ Joins lift according to rules from paper

        lift: output of lempels_lift()
        digit: number of digits added to longest strings of same character
            before lifting; this is relevant for determining how to join

        Return value: lifted sequence as a list of integers
    """
    if len(lift) == 1:
        # Just one string from lift; easiest case
        seq = lift[0]
    elif digit:
        # Easy case where the join is performed on a lift of the longest string
        # of 1s
        seq = easy_join(lift, ceil(log(total_length, K)))
    else:
        # Hard case where the join is performed on lifts of the longest string
        # of 1s plus an unspecified character on either side
        assert len(lift) == K # This is true by a lemma in the paper

    return seq

def pkl_via_lempels_lift(K, L):
    """ Generates a P(K)L-sequence using Lempel's lift.

        K: alphabet size
        L: desired length

        Return value: P(K)L-sequence as a list of integers
    """
    # Decompose into base-K representation
    digits = length_in_base(L, K)
    seq = deque(range(1, digit + 1))
    for i, digit in enumerate(digits[1:]):
        # Lift and join
        seq = joined_lift(lempels_lift(seq, K))
        # Lengthen seq
        if digit:
            seq_length = len(seq)
            streak_length = floor(log(seq_length, K))
            streak = 0
            for j in range(seq_length - streak_length + 2, seq_length):
                if seq[j] % K == seq[j-1] % K:
                    streak += 1
                else:
                    streak = 0
            j = 0
            while j < seq_length + digit:
                if seq[j] % K == seq[j-1] % K:
                    streak += 1
                else:
                    streak = 0
                if streak == streak_length and 1 <= seq[j] <= digit:
                    seq.insert(seq[j], j + 1)
                    j += 1
                j += 1
    return seq

def help_formatter(prog):
    """ So formatter_class's max_help_position can be changed. """
    return argparse.HelpFormatter(prog, max_help_position=40)

if __name__ == "__main__":
    import argparse
    from sys import stdin, stderr
    _intro = ("Specify either 'check' or 'construct' to check that a sequence "
              "is a P(K)L-sequence or construct one, respectively.")
    parser = argparse.ArgumentParser(
            description=_intro,
            formatter_class=help_formatter
        )
    subparsers = parser.add_subparsers(
            help=("construct a P(K)L-sequence or check that a sequence is a "
                  "P(K)L-sequence"),
            dest="subparser_name"
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
    if args.subparser_name is None:
        print(_intro, file=stderr)
    if args.subparser_name == "check":
        with stdin as input_stream:
            if not args.separator:
                seq = list(input_stream.read().strip())
            else:
                input_stream.read().strip().split(args.separator)
        if not check(seq, low_mem=args.low_memory, temp_dir=args.temp_dir):
            raise RuntimeError('Input is not a P(K)L-sequence.')
        print("Input is a P(K)L-sequence.", file=stderr)
    else:
        assert args.subparser_name == "construct"
