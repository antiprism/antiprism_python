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
Make a model with axial symmetry based on a belt of staggered pentagons
'''

import argparse
import sys
import math
import anti_lib
from anti_lib import Vec

epsilon = 1e-13
phi = (math.sqrt(5) + 1)/2


def calculate_belt_points(pgon, model_type):
    ang = pgon.angle()/2
    tan_a = math.sin(ang)/(math.cos(ang)+phi)     # half of unit edge
    P = Vec(1/(2*tan_a), -1/2, 0)
    tan_b = math.sin(ang)/(math.cos(ang)+1/phi)   # half of phi edge
    Q = Vec(phi/(2*tan_b), -phi/2, 0)
    bar_ht = math.cos(math.pi/10)           # planar height of pentagon "bar"
    diff_r = P[0] - Q[0]
    try:
        ht0 = math.sqrt(bar_ht**2 - diff_r**2)/2
    except:
        raise ValueError('model is not constructible')

    P[2] = -ht0
    Q[2] = ht0
    pent_ht = math.sqrt(5 + 2*math.sqrt(5))/2   # planar height of pentagon
    x_cap = P[0] - pent_ht*(diff_r/bar_ht)
    z_cap = -ht0 + pent_ht*((2*ht0)/bar_ht)
    R = Vec(x_cap, 0, z_cap)
    S = Vec(x_cap*math.cos(ang), x_cap*math.sin(ang), -z_cap)
    cap_inrad = x_cap*math.cos(ang)
    apex_ht = 0
    if model_type == 'a':
        try:
            apex_ht = z_cap + cap_inrad * (z_cap + P[2])/(P[0] - cap_inrad)
        except:
            raise ValueError('could not calculate apex height')
    A = Vec(0, 0, apex_ht)
    return [P, Q, R, S, A]


def make_model(pgon, model_type):
    pts = calculate_belt_points(pgon, model_type)
    N = pgon.N
    ang = pgon.angle()
    points = []
    for point in pts[0:2]:
        for i in range(N):
            points.append(point.rot_z(i*ang))
            points.append((point - Vec(0, 2*point[1], 0)).rot_z(i*ang))

    for point in pts[2:4]:
        points += [point.rot_z(i*ang) for i in range(N)]

    if model_type == 'a':
        points += [pts[4], -pts[4]]

    faces = []
    if model_type == 't':
        for off in [4*N, 5*N]:
            faces += [[i+off for i in range(N)]]

    for i in range(N):
        faces += [[2*N + 2*i, 4*N + i, 2*N + 2*i+1, 2*i+1, 2*i]]
        faces += [[2*i+1, 5*N + i, (2*i+2) % (2*N),
                   2*N + (2*i+2) % (2*N), 2*N + 2*i+1]]
        if model_type == 't':
            faces += [[2*i, 2*i+1, 5*N + i, 5*N + (N+i-1) % N]]
            faces += [[2*N + 2*i+1, 2*N + (2*i+2) % (2*N),
                       4*N + (i+1) % N, 4*N + i]]
        else:
            faces += [[2*i, 2*i+1, 5*N + i, 6*N + 1, 5*N + (N+i-1) % N]]
            faces += [[2*N + 2*i+1, 2*N + (2*i+2) % (2*N),
                       4*N + (i+1) % N, 6*N, 4*N + i]]

    return points, faces


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'polygon',
        help='number of pairs of pentagons (default: 4) '
             '(or may be a polygon fraction, e.g. 5/2)',
        type=anti_lib.read_polygon,
        nargs='?',
        default='7')
    parser.add_argument(
        '-m', '--model_type',
        help='model type: a - with apex, t - truncated apex',
        default='t')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    pgon = args.polygon

    try:
        points, faces = make_model(pgon, args.model_type)
    except Exception as e:
        parser.error(e.args[0])

    out = anti_lib.OffFile(args.outfile)
    out.print_all_pgon(points, faces, pgon)

if __name__ == "__main__":
    main()
