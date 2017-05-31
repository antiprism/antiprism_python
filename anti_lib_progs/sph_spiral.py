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
Distribute points in a spiral on a sphere.
'''

import argparse
import sys
import math
import anti_lib
from anti_lib import Vec


def angle_to_point(a, number_turns):
    a2 = 2 * a * number_turns   # angle turned around y axis
    r = math.sin(a)             # distance from y axis
    y = math.cos(a)
    x = r * math.sin(a2)
    z = r * math.cos(a2)
    return Vec(x, y, z)


# binary search for angle with distance rad from point with angle a0
def psearch(a1_delt, a1, a0, rad, number_turns):
    a_test = a1_delt + (a1 - a1_delt) / 2.0
    dist = (angle_to_point(a_test, number_turns) -
            angle_to_point(a0, number_turns)).mag()
    eps = 1e-5
    if rad + eps > dist > rad - eps:
        return a_test
    elif rad < dist:      # Search in first interval
        return psearch(a1_delt, a_test, a0, rad, number_turns)
    else:                 # Search in second interval
        return psearch(a_test, a1, a0, rad, number_turns)


def calc_points(args):
    points = []
    number_turns = args.number_turns
    if not number_turns:
        number_turns = 1e-12
    if args.distance_between_points:
        rad = 2*args.distance_between_points
    else:
        # half distance between turns on a rad 1 sphere
        rad = 2*math.sqrt(1 - math.cos(math.pi/(number_turns-1)))
    a0 = 0
    cur_point = Vec(0.0, 1.0, 0.0)
    points.append(cur_point)

    delt = math.atan(rad / 2) / 10
    a1 = a0 + .0999999 * delt                        # still within sphere
    while a1 < math.pi:
        if (cur_point - angle_to_point(a1, number_turns)).mag() > rad:
            a0 = psearch(a1 - delt, a1, a0, rad, number_turns)
            cur_point = angle_to_point(a0, number_turns)
            points.append(cur_point)
            a1 = a0

        a1 += delt

    return points


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'number_turns',
        help='number of times the spiral turns around the '
             'axis (default: 10)',
        type=float,
        nargs='?',
        default=10)
    parser.add_argument(
        'distance_between_points',
        help='the distance between consequetive points '
             '(default: the distance beteen spiral turns)',
        type=anti_lib.read_positive_float,
        nargs='?')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    points = calc_points(args)

    out = anti_lib.OffFile(args.outfile)
    out.print_all(points, [])

if __name__ == "__main__":
    main()
