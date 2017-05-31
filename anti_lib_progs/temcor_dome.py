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
Make a Temcor-style dome, using the method described in
https://groups.google.com/d/msg/geodesichelp/hJ3V9Nfp3kE/nikgoBPSFfwJ .
The base model is a pyramid with a unit edge base polygon
at a specified height above the origin. The axis to rotate the plane
about passes through the origin and is in the direction of the base
polygon mid-edge to the pyramid apex.
'''

import argparse
import sys
import anti_lib
from anti_lib import Vec, Mat


def calc_temcor_side(pgon, pyramid_ht, base_ht, freq):
    freq += 1
    inrad = pgon.inradius()
    axis = Vec(-inrad, 0, pyramid_ht)           # axis to rotate plane around
    A = Vec(0, 0, base_ht + pyramid_ht)         # apex
    B = Vec(inrad, 0.5, base_ht)                # base polygon vertex

    n0 = Vec.cross(A, B).unit()

    edge_ang = anti_lib.angle_around_axis(Vec(B[0], 0, B[2]), B, axis)
    ang_inc = edge_ang/freq

    points = []
    faces = []
    for i in range(freq):
        n_inc = Mat.rot_axis_ang(axis, i*ang_inc) * Vec(0, 1, 0)
        edge_v = Vec.cross(n_inc, n0).unit()
        last_idx = i*(i-1)//2
        new_idx = i*(i+1)//2
        for j in range(i + 1):
            v = Mat.rot_axis_ang(axis, -2*j*ang_inc) * edge_v
            points.append(v)
            if new_idx and j < i:
                faces.append([new_idx+j, new_idx+j+1, last_idx+j])
            if j < i-1:
                faces.append([new_idx+j+1, last_idx+j+1, last_idx+j])

    return (points, faces)


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'number_sides',
        help='number of sides (default: 6) (or may be a '
             'polygon fraction, e.g. 5/2)',
        type=anti_lib.read_polygon,
        nargs='?',
        default='6')
    parser.add_argument(
        'pyramid_ht',
        help='pyramid_height (default: 0.5)',
        type=float,
        nargs='?',
        default=0.5)
    parser.add_argument(
        'base_ht',
        help='base_height (default: 0.5)',
        type=float,
        nargs='?',
        default=0.5)
    parser.add_argument(
        'frequency',
        help='frequency (default: 5)',
        type=anti_lib.read_positive_int,
        nargs='?',
        default=5)
    parser.add_argument(
        '-t', '--triangle-only',
        help='only output one triangle section of the dome',
        action='store_true')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    if abs(args.pyramid_ht + args.base_ht) < 1e-13:
        parser.error('base height cannot be the negative of the '
                     'pyramid height')

    (points, faces) = calc_temcor_side(args.number_sides, args.pyramid_ht,
                                       args.base_ht, args.frequency)

    out = anti_lib.OffFile(args.outfile)
    if args.triangle_only:
        out.print_all(points, faces)
    else:
        out.print_all_pgon(points, faces, args.number_sides,
                           repeat_side=True)

if __name__ == "__main__":
    main()
