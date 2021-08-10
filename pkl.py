#!/usr/bin/env python3
"""
pkl.py

Tools for constructing P(K)L-sequences as described in

A. Nellore and R. Ward. "Arbitrary-length analogs to de Bruijn sequences"

and verifying that constructed sequences are P(K)L-sequences.

The implementation here is didactic. While it is uses O(L log K) space and
takes O(L) time, as discussed in the paper above, several optimizations could
be made.

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
from math import ceil, floor, gcd
from math import log as mathlog
from os import remove
from collections import deque, defaultdict
from itertools import combinations

def log(x, b):
    """ Gets around Python's precision issues when computing logs.

        x: integer
        b: integer

        Return value: log_x b
    """
    y = b
    power = 1
    while y < x:
        y *= b
        power += 1
    if y == x:
        return power
    return mathlog(x, b)

def check(seq, low_mem=False, temp_dir=None):
    """ Checks if a (circular) string is a P(K)L-sequence.

        seq: string or list of integers to check
        low_mem: True iff check should use sqlite dictionary on disk rather
            than memory
        temp: path to temporary directory or None if default should be used

        Return value: True iff seq is a P(K)L-sequence.
    """
    L = len(seq)
    K = max(len(set(seq)), 2)
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
    for j in range(K // repetitions):
        lift.append([0 for _ in range(repetitions * L)])
        lift[-1][-1] = j
        for i in range(L * repetitions):
            lift[-1][i] = (lift[-1][i-1] + seq[i % L]) % K
    return lift

def easy_hamiltonian_path(lift, streak_length, K):
    """ Writes path for joining lift on (x, x + 1, ..., x + streak_length - 1)

        streak_length: length of join string
        lift: strings to join, where the nth string is 1 + the (n-1)th string
        K: alphabet size

        Return value: deque of tuples, each representing an edge of a path and
            of the form

            (index of start vertex, index of end vertex,
             start pos of join substring in start vertex string,
             start pos of join substring in end vertex string).

            Edges are in order of their appearance in the path.
    """
    string_length = len(lift[0])
    total_strings = len(lift)
    streak = 1
    # So streaks can be detected starting at first character
    for i in range(string_length - streak_length + 2, string_length):
        if lift[0][i % string_length] % K == (lift[0][
                                                    (i - 1) % string_length
                                                ] + 1) % K:
            streak += 1
        else:
            streak = 1
    for i in range(string_length):
        if lift[0][i % string_length] % K == (lift[0][
                                                    (i - 1) % string_length
                                                ] + 1) % K:
            streak += 1
        else:
            streak = 1
        if streak == streak_length:
            return deque([
                    (j, j + 1, (i - streak + 2) % string_length,
                        (i - streak + 1) % string_length)
                        for j in range(total_strings - 1)
                ])

def cycle_join(strings, spanning_trees):
    """ Joins circular strings in a Lempel's lift according to spanning trees

        Spanning trees are of "join graph" from function description of
        spanning_trees(). SIDE EFFECT: spanning_trees gets eaten

        strings: list of deques representing circular strings to join
        spanning_trees: deque of tuples, each representing an edge and of the
            form

            (index of start vertex, index of end vertex,
             start pos of join substring in start vertex string,
             start pos of join substring in end vertex string).

            Edges should be ordered as the result of appending depth-first
            traversals of the spanning trees to a deque.

        Return value: list of lists, each giving a joined string for a
            spanning tree
    """
    # For keeping track of where joins of previously added strings are
    substring_displacements = defaultdict(lambda: defaultdict(int))
    # Initialize all displacements
    for edge in spanning_trees:
        substring_displacements[edge[0]][edge[2]] = edge[2]
        substring_displacements[edge[1]][edge[3]] = edge[3]
    joined = []
    last_head = None
    seen = set()
    while spanning_trees:
        tail, head, tail_start, head_start = spanning_trees.popleft()
        if tail not in seen:
            # Start new string
            joined.append(deque(strings[tail][tail_start:]
                          + strings[tail][:tail_start]
                          + strings[head][head_start:]
                          + strings[head][:head_start]))
            tail_string_length = len(strings[tail])
            for start in substring_displacements[tail]:
                substring_displacements[tail][start] = (
                        substring_displacements[tail][start] - tail_start
                    ) % tail_string_length
            for start in substring_displacements[head]:
                substring_displacements[head][start] = tail_string_length + (
                        substring_displacements[head][start] - head_start
                    ) % len(strings[head])
            seen.update([tail, head])
        else:
            # Rotate string and concatenate
            shift = -substring_displacements[tail][tail_start]
            joined[-1].rotate(shift)
            length = len(joined[-1])
            for vertex in seen:
                for pos in substring_displacements[vertex]:
                    substring_displacements[vertex][pos] = (
                            substring_displacements[vertex][pos] + shift
                        ) % length
            joined[-1].extend(strings[head][head_start:])
            joined[-1].extend(strings[head][:head_start])
            seen.add(head)
            for start in substring_displacements[head]:
                substring_displacements[head][start] = length + (
                        substring_displacements[head][start] - head_start
                    ) % len(strings[head])
    return list(map(list, joined))

def hard_spanning_trees(strings, streak_length, K, with_constant=False):
    """ Gives spanning trees of connected components of join graph of strings

        The "join graph" of strings is a graph where each vertex represents a
        string from strings, and each (undirected) edge is labeled by

        (when with_constant=True)
        a substring of the form (y, x, x + 1, ..., x + streak_length - 1) or
        (x, x + 1, ..., x + streak_length - 1, y)

        (when with_constant=False)
        a substring of the form (x, x+1, ..., x + streak_length - 1)

        at which a join can be performed between the strings corresponding to
        the vertices it connects.

        strings: list of deques representing circular strings
        streak_length: length of streak in substring
            (x, x + 1, ... , x + streak_length - 1)
        K: size of alphabet

        Return value: deque of tuples, each representing an edge of a spanning
            tree and of the form

            (index of start vertex, index of end vertex,
             start pos of join substring in start vertex string,
             start pos of join substring in end vertex string).

            Edges are ordered from appending depth-first traversals of the
            spanning trees to a deque.
    """
    string_count = len(strings)
    if with_constant:
        seen_next = defaultdict(list) # dictionary mapping first character of
                                      # streak and character seen AFTER streak
                                      # to list of tuples, each in format
                                      # (vertex index, substring position)
        seen_prev = defaultdict(list) # dictionary mapping first character of
                                      # streak and character seen BEFORE streak
                                      # to list of tuples, each in format
                                      # (vertex index, substring position)
        for i in range(string_count):
            string_length = len(strings[i])
            streak = 1
            # So streaks can be detected starting at first character
            for j in range(string_length - streak_length + 2, string_length):
                if strings[i][j % string_length] % K == (strings[i][
                        (j - 1) % string_length
                    ] + 1) % K:
                    streak += 1
                else:
                    streak = 1
            for j in range(string_length):
                if strings[i][j % string_length] % K == (strings[i][
                        (j - 1) % string_length
                    ] + 1) % K:
                    streak += 1
                else:
                    streak = 1
                if streak == streak_length:
                    seen_next[
                            (strings[i][(j - streak + 1) % string_length],
                             strings[i][(j + 1) % string_length])
                        ].append((i, (j - streak + 1) % string_length))
                    seen_prev[
                            (strings[i][(j - streak + 1) % string_length],
                             strings[i][(j - streak) % string_length])
                        ].append((i, (j - streak) % string_length))
                elif streak > streak_length:
                    assert streak_length == streak - 1
                    seen_next[
                            (strings[i][(j - streak + 2) % string_length],
                             strings[i][(j + 1) % string_length])
                        ].append((i, (j - streak + 2) % string_length))
                    seen_prev[
                            (strings[i][(j - streak + 2) % string_length],
                             strings[i][(j - streak + 1) % string_length])
                        ].append((i, (j - streak + 1) % string_length))
        # Construct graph
        graph = [[] for _ in range(string_count)] # List of lists of 
                                                  # (head index,
                                                  # start pos of tail,
                                                  # start pos of head)
        for seen in (seen_prev, seen_next):
            for substring in seen:
                for combo in combinations(seen[substring], 2):
                    graph[combo[0][0]].append(
                            (combo[1][0], combo[0][1], combo[1][1])
                        )
    else:
        seen = defaultdict(list) # dictionary mapping first character of streak
                                 # to list of tuples, each in the format
                                 # (vertex index, substring position)
        for i in range(string_count):
            string_length = len(strings[i])
            streak = 1
            # So streaks can be detected starting at first character
            for j in range(string_length - streak_length + 2, string_length):
                if strings[i][j % string_length] % K == (strings[i][
                        (j - 1) % string_length
                    ] + 1) % K:
                    streak += 1
                else:
                    streak = 1
            for j in range(string_length):
                if strings[i][j % string_length] % K == (strings[i][
                        (j - 1) % string_length
                    ] + 1) % K:
                    streak += 1
                else:
                    streak = 1
                if streak == streak_length:
                    seen[strings[i][(j - streak + 1) % string_length]].append(
                            (i, (j - streak + 1) % string_length)
                        )
                elif streak > streak_length:
                    assert streak_length == streak - 1
                    seen[strings[i][(j - streak + 2) % string_length]].append(
                            (i, (j - streak + 2) % string_length)
                        )
        # Construct graph
        graph = [[] for _ in range(string_count)] # List of lists of
                                                  # (head index,
                                                  # start pos of tail,
                                                  # start pos of head)
        for substring in seen:
            for combo in combinations(seen[substring], 2):
                graph[combo[0][0]].append(
                        (combo[1][0], combo[0][1], combo[1][1])
                    )
    # Obtain spanning trees via depth-first search
    not_seen = set(range(string_count))
    stack = deque()
    traversal = deque()
    while True:
        try:
            tail = stack.pop()
        except IndexError:
            try:
                tail = not_seen.pop()
            except KeyError:
                # Everything's been seen
                break
        for head in graph[tail]:
            if head[0] in not_seen:
                stack.append(head[0])
                not_seen.remove(head[0])
                traversal.append((tail,) + head)
    return traversal

def joined_lift(lift, digit, K):
    """ Joins lift according to rules from paper

        lift: output of lempels_lift()
        digit: number of digits added to longest strings of same character
            before lifting; this is relevant for determining how to join
        K: alphabet size

        Return value: lifted sequence as a list of integers
    """
    if len(lift) == 1:
        # Just one string from lift; easiest case
        seq = lift[0]
    else:
        total_length = sum(map(len, lift))
        floorlog = floor(log(total_length, K))
        ceillog = ceil(log(total_length, K))
        if digit or floorlog == ceillog:
            # Easy case where the join is performed on a lift of the longest
            # string of 1s
            seq = cycle_join(lift, easy_hamiltonian_path(lift, ceillog, K))[0]
        else:
            # Hard case where the join is performed on lifts of the longest
            # string of 1s plus an unspecified character on either side
            trees = hard_spanning_trees(
                        lift, floorlog - 1, K, with_constant=True
                    )
            if trees:
                joined = cycle_join(lift, trees)
                if len(joined) > 1:
                    # Need to merge on string of form
                    # (x, x + 1, ... x + streak_length - 1)
                    joined = cycle_join(joined, hard_spanning_trees(
                            joined, floorlog - 1, K, with_constant=False
                        ))
                    assert len(joined) == 1
                else:
                    assert joined
            else:
                joined = cycle_join(lift, hard_spanning_trees(
                            lift, floorlog - 1, K, with_constant=False
                        ))
                assert len(joined) == 1
            seq = joined[0]
    return seq

def pkl_via_lempels_lift(K, L):
    """ Generates a P(K)L-sequence using Lempel's lift.

        K: alphabet size
        L: desired length

        Return value: P(K)L-sequence as a deque of integers
    """
    # Decompose into base-K representation
    digits = length_in_base(L, K)
    seq = list(range(1, digits[0] + 1))
    for i, digit in enumerate(digits[1:]):
        # Lift and join
        seq = joined_lift(lempels_lift(seq, K), digits[i], K)
        # Lengthen seq
        if digit:
            seq_length = len(seq)
            streak_length = floor(log(seq_length, K))
            new_seq = []
            streak = 1
            for j in range(seq_length - streak_length + 2, seq_length):
                if seq[j] % K == seq[j-1] % K:
                    streak += 1
                else:
                    streak = 1
            dones = set()
            for j in range(seq_length):
                if seq[j % seq_length] % K == seq[(j-1) % seq_length] % K:
                    streak += 1
                else:
                    streak = 1
                new_seq.append(seq[j])
                if (streak == streak_length
                    and 1 <= seq[j] <= digit and seq[j] not in dones):
                    new_seq.append(seq[j])
                    dones.add(seq[j])
                    if len(dones) == digit:
                        seq = new_seq + seq[j+1:]
                        break
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
        quit()
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
        seq = pkl_via_lempels_lift(args.alphabet_size, args.sequence_length)
        print(seq)
        if check(seq):
            print("Output is a P(K)L-sequence.", file=stderr)
