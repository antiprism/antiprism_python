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
Create a polyhedron which tiles the sphere with congruent triangles.
'''

import argparse
import sys
import math
import anti_lib
from anti_lib import Vec, Mat


def get_base_points(pgon, centres=0, twist=False):
    N = pgon.N
    D = pgon.D
    a = pgon.angle()/2
    points = [Vec(0, 0, 1), Vec(0, 0, -1)]  # poles
    A = a
    B = math.pi/2 * (1 - D/N)
    cos_lat = 1/(math.tan(A) * math.tan(B))
    sin_lat = math.sqrt(1 - cos_lat**2)
    p_north = Vec(sin_lat, 0, cos_lat)
    p_south = Vec(sin_lat, 0, -cos_lat)
    for i in range(N):
        points.append(p_north.rot_z(i*2*a))
    for i in range(N):
        points.append(p_south.rot_z(i*2*a+a))
    if centres:
        cos_cent_lat = math.cos(B)/math.sin(A)
        sin_cent_lat = math.sqrt(1 - cos_cent_lat**2)
        cent_north = Vec(sin_cent_lat, 0, cos_cent_lat)
        cent_south = Vec(sin_cent_lat, 0, -cos_cent_lat)
        for i in range(N):
            points.append(cent_north.rot_z(i*2*a+a) * centres)
        for i in range(N):
            points.append(cent_south.rot_z(i*2*a) * centres)

    rhombi = []
    for i in range(N):
        rhombi.append([2 + i, 0, 2 + ((i+1) % N), 2 + N + i])
    for i in range(N):
        rhombi.append([2 + N + i, 1, 2 + N + ((i-1) % N), 2 + i])

    if twist:
        hx = [0, 2, 2 + N, 1, 2 + N+N//2, 2 + N-N//2]
        axis = points[hx[0]] + points[hx[2]] + points[hx[4]]
        rot = Mat.rot_axis_ang(axis.unit(), 2*math.pi/3)
        for i in list(range(hx[5]+1, 2 + N)) + list(range(hx[4]+1, 2 + 2*N)):
            points[i] = rot * points[i]
        if centres:
            for r in list(range(N-N//2, N+1)) + list(range(2*N-N//2, 2*N)):
                points[2 + 2*N + r] = rot * points[2 + 2*N + r]
        for r in list(range(N-N//2, N+1)) + list(range(2*N-N//2, 2*N)):
            for i in range(4):
                for idx, p_idx in enumerate(hx):
                    if p_idx == rhombi[r][i]:
                        rhombi[r][i] = hx[(idx+2) % 6]
                        break

    return points, rhombi


def make_gyroelongated_dipyramid_model(rhombi):
    faces = []
    for r in rhombi:
        faces.append([r[0], r[1], r[2]])
        faces.append([r[2], r[3], r[0]])
    return faces


def make_scalahedron_model(rhombi):
    faces = []
    for r in rhombi:
        faces.append([r[1], r[2], r[3]])
        faces.append([r[3], r[0], r[1]])
    return faces


def make_subscalahedron_model(rhombi):
    faces = []
    for idx, r in enumerate(rhombi):
        apex = 2 + len(rhombi) + idx
        for i in range(4):
            faces.append([r[i], r[(i+1) % 4], apex])
    return faces


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'polygon_fraction',
        help='number of sides of the base polygon (N), '
        'or a fraction for star polygons (N/D) (default: 5)',
        default='5',
        nargs='?',
        type=anti_lib.read_polygon)
    parser.add_argument(
        'poly_type',
        help='polyhedron type: '
             '0 - base rhombic tiling, '
             '1 - gyroelongated dipyramid, '
             '2 - scalenohedron, '
             '3 - subscalenohedron '
             '4 - inverted subscalenohedron (default: 1)',
        choices=['0', '1', '2', '3', '4'],
        default='1')
    parser.add_argument(
        '-t', '--twist',
        help='twist one half of the model by 120 degrees'
             '(N must be odd)',
        action='store_true')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()
    pgon = args.polygon_fraction
    if args.twist:
        msg_start = 'twist: cannot twist model, polygon fraction '
        if pgon.parts > 1:
            parser.error(msg_start + 'is compound')
        msg = ''
        if not pgon.N % 2:
            msg = 'numerator'
        elif not pgon.D % 2:
            msg = 'denominator'
        if msg:
            parser.error(msg_start + msg + ' is even')

    if args.poly_type == '3':
        add_centres = 1
    elif args.poly_type == '4':
        add_centres = -1
    else:
        add_centres = 0

    points, rhombi = get_base_points(pgon, add_centres, args.twist)
    if args.poly_type == '0':
        faces = rhombi
    elif args.poly_type == '1':
        faces = make_gyroelongated_dipyramid_model(rhombi)
    elif args.poly_type == '2':
        faces = make_scalahedron_model(rhombi)
    elif args.poly_type == '3' or args.poly_type == '4':
        faces = make_subscalahedron_model(rhombi)
    else:
        parser.error('unknown polyhedron type')

    out = anti_lib.OffFile(args.outfile)
    out.print_all_pgon(points, faces, pgon)

if __name__ == "__main__":
    main()
