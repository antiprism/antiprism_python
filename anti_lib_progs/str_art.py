#!/usr/bin/env python3

# Copyright (c) 2014-2016 Adrian Rossiter <adrian@antiprism.com>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

'''
Simple epicycloid string art patterns. Output is in OFF format.
'''

import argparse
import sys
import math
import fractions


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'factor',
        help='an integer or fraction greater than 1',
        type=fractions.Fraction,
        nargs='?',
        default=2)
    parser.add_argument(
        'num_pins',
        help='number of pins',
        type=int,
        nargs='?',
        default=120)
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()
    if args.factor < 1:
        parser.error('factor cannot be less than 1')
    if args.num_pins < 2:
        parser.error('number of pins cannot be less than 2')

    N = args.factor.numerator
    D = args.factor.denominator
    pins = args.num_pins
    turns = N / fractions.gcd(N, D)
    strings = turns*pins
    print('OFF\n{} {} 0'.format(pins, strings), file=args.outfile)

    for i in range(pins):
        print(math.cos(2*math.pi*i/pins), math.sin(2*math.pi*i/pins), 0,
              file=args.outfile)

    for i in range(pins):
        print(2, (i * N) % pins, (i * D) % pins, file=args.outfile)

if __name__ == "__main__":
    main()
