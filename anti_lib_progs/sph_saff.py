#!/usr/bin/env python3

# Copyright (c) 2003-2016 Adrian Rossiter <adrian@antiprism.com>
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
Distribute num_points (default 20) on a sphere using the algorithm from
"Distributing many points on a sphere" by E.B. Saff and
A.B.J. Kuijlaars, Mathematical Intelligencer 19.1 (1997) 5--11.
'''

import argparse
import sys
import math
import anti_lib
from anti_lib import Vec


def calc_points(args):
    points = []
    use_angle = 'angle' in args and args.angle is not None
    if use_angle:
        ang = (args.angle * math.pi/180) % (2*math.pi)

    N = args.number_points
    for k in range(1, N + 1):
        h = -1 + 2 * (k - 1) / float(N - 1)
        theta = math.acos(h)
        if k == 1 or k == N:
            phi = 0
        elif use_angle:
            phi += ang
        else:
            phi += 3.6 / math.sqrt(N * (1 - h * h))

        points.append(Vec(math.sin(phi) * math.sin(theta),
                      math.cos(phi) * math.sin(theta),
                      -math.cos(theta)))
        phi %= 2*math.pi

    return points


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'number_points',
        help='number of points to distribute on a sphere',
        type=anti_lib.read_positive_int,
        nargs='?',
        default=100)
    parser.add_argument(
        '-a', '--angle',
        help='increment each point placement by a fixed angle instead '
             'of using the Saff and Kuiljaars placement method',
        type=float)
    parser.add_argument(
        '-x', '--exclude-poles',
        help='exclude the pole point circles',
        action='store_true')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    points = calc_points(args)
    start = int(args.exclude_poles)
    end = len(points)-start

    out = anti_lib.OffFile(args.outfile)
    out.print_all(points[start:end], [])

if __name__ == "__main__":
    main()
