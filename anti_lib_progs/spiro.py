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
Spirograph generator. Output is in OFF format.
'''

import argparse
import sys
import math
import fractions


def spiro(num_teeth_fixed, num_teeth_move, height, num_segs, outfile):
    N = abs(num_teeth_fixed)
    D = abs(num_teeth_move)
    side_sign = num_teeth_move/D
    height = height*D/N
    num_segs = num_segs

    turns = D/fractions.gcd(N, D)
    print('OFF\n{} 1 0'.format(num_segs), file=outfile)

    for i in range(num_segs):
        ang_fixed = 2*math.pi*turns*i/num_segs
        ang_move = side_sign * ang_fixed * N/D
        move_cent = [math.cos(ang_fixed)*(N + side_sign*D)/N,
                     math.sin(ang_fixed)*(N + side_sign*D)/N, 0]
        move_offset = [height*math.cos(ang_fixed+ang_move),
                       height*math.sin(ang_fixed+ang_move), 0]
        P = [move_cent[i] + move_offset[i] for i in range(3)]
        print(P[0], P[1], P[2], file=outfile)

    print(num_segs, *[i for i in range(num_segs)], file=outfile)


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'num_teeth_fixed',
        help='number of teeth of fixed gear wheel',
        type=int,
        nargs='?',
        default=20)
    parser.add_argument(
        'num_teeth_move',
        help='number of teeth on moving gear wheel (positive to run on '
             'outside fixed wheel or negative to run on inside)',
        type=int,
        nargs='?',
        default=-8)
    parser.add_argument(
        'height',
        help='height of drawing point on moving gear from its centre, '
             'as ratio of moving gear radius',
        type=float,
        nargs='?',
        default=1.0)
    parser.add_argument(
        '-n', '--num-segs',
        help='number of segments used to draw pattern',
        type=int,
        default=1000)
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    if args.num_teeth_fixed < 2:
        parser.error('number of teeth on fixed gear wheel cannot be '
                     'less than 2')

    if abs(args.num_teeth_move) < 2:
        parser.error('number of teeth on moving gear wheel cannot be '
                     'less than 2')

    if args.num_segs < 2:
        parser.error('number of segments cannot be less than 2')

    spiro(args.num_teeth_fixed, args.num_teeth_move,
          args.height, args.num_segs, args.outfile)

if __name__ == "__main__":
    main()
