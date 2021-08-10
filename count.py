#!/usr/bin/env python3
"""
count.py

Counts P(K)L sequences for any value of K and L. Generates necklaces
very inefficiently; an improvement would implement, e.g.,

    A fast algorithm to generate necklaces with fixed content,
    Theoretical Computer Science archive Volume 301, Issue 1-3 (May 2003)
    doi:10.1016/S0304-3975(03)00049-5

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

import argparse
import math
import itertools
from collections import defaultdict
from pkl import check

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequence-length', '-l', required=True, type=int,
            help='length of P(K)L sequence'
        )
    parser.add_argument('--alphabet-size', '-k', required=True, type=int,
            help='alphabet size'
        )
    args = parser.parse_args()
    count = 0
    for seq in itertools.product(
                *([list(range(args.alphabet_size))]*args.sequence_length)
            ):
        if check(seq):
            count += 1
    print("{}\t{}\t{}".format(
            args.alphabet_size,
            args.sequence_length,
            count // args.sequence_length)
        )
