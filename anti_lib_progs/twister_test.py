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
Twist two polygons placed on axes at a specified angle and joined by a vertex
'''


import argparse
import sys
import math
import anti_lib
from anti_lib import Vec, Mat

epsilon = 1e-13


def calc_polygons2(pgon0, r0, ang, pgon1, r1, angle_between_axes):
    # rotate initial vertex of polygon0 by ang. Common
    # point lies on a line through vertex in the direction
    # of axis0 (z-axis). Common vertex lies on a cylinder
    # radius r1 with axis on axis1.

    # rotate model so axis1 is on z-axis, to simplify the
    # calculation of the intersection
    rot = Mat.rot_axis_ang(Vec(0, 1, 0), angle_between_axes)

    # initial vertex turned on axis
    V = Vec(r0, 0, 0).rot_z(ang)
    # vertex in rotated model
    Q = rot * V
    # direction of line in rotated model
    u = rot * Vec(0, 0, 1)

    # equation of line is P = Q + tu for components x y z and parameter t
    # equation of cylinder at z=0 is x^2 + y^2 = r1^2
    a = u[0]**2 + u[1]**2
    b = 2*(Q[0]*u[0] + Q[1]*u[1])
    c = Q[0]**2 + Q[1]**2 - r1**2

    disc = b**2 - 4*a*c
    if disc < -epsilon:
        raise Exception("model is not geometrically constructible")
    elif disc < 0:
        disc = 0

    t = (-b - math.sqrt(disc))/(2*a)   # negative gives most distant point

    P = V + Vec(0, 0, t)   # the common point

    points = []
    points += [P.rot_z(i*pgon0.angle()) for i in range(pgon0.N)]
    faces = [[i for i in range(pgon0.N)]]

    Q = rot * P
    rot_inv = Mat.rot_axis_ang(Vec(0, 1, 0), -angle_between_axes)
    points += [rot_inv * Q.rot_z(i*pgon1.angle()) for i in range(pgon1.N)]
    faces += [[i+pgon0.N for i in range(pgon1.N)]]

    return points, faces


def read_turn_angle(ang_str):
    if ang_str[-1] == 'e':
        ang_type = 'e'
        ang_str = ang_str[:-1]
    elif ang_str[-1] == 'x':
        ang_type = 'x'
        ang_str = ang_str[:-1]
    else:
        ang_type = 'd'
    try:
        ang_val = float(ang_str)
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            'invalid numeric value: ', e.args[0])

    return [ang_val, ang_type]


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'number_sides0',
        help='number of sides of polygon 0 (default: 6) '
             '(or may be a polygon fraction, e.g. 5/2)',
        type=anti_lib.read_polygon,
        nargs='?',
        default='6')
    parser.add_argument(
        'number_sides1',
        help='number of sides of polygon 0 (default: 5) '
             '(or may be a polygon fraction, e.g. 5/2)',
        type=anti_lib.read_polygon,
        nargs='?',
        default='5')
    parser.add_argument(
        '-A', '--angle_between_axes',
        help='angle between the two axes (default: 60 degs)',
        type=float,
        default=60.0)
    parser.add_argument(
        '-r', '--ratio',
        help='ratio of edge lengths (default: 1.0)',
        type=float,
        default=1.0)
    parser.add_argument(
        '-a', '--angle',
        help='amount to turn polygon on axis0 in degrees '
        '(default: 0), or a value followed by \'e\', '
             'where 1.0e is half the central angle of an edge, '
             'which produces an edge-connected model (negative '
             'values may have to be specified as, e.g. -a=-1.0e), '
             'or a value followed by x, which is like e but with '
             'a half turn offset',
        type=read_turn_angle,
        default='0')
    parser.add_argument(
        '-x', '--x-axis-vert',
        help='offset of vertex of side polygon to align with x-axis, with 0'
             'being the vertex attached to the axial polygon',
        type=int)
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    pgon0 = args.number_sides0
    pgon1 = args.number_sides1

    if args.angle[1] == 'e':     # units: half edge central angle
        turn_angle = args.angle[0] * pgon0.angle()/2
    elif args.angle[1] == 'x':   # units: half edge central angle
        turn_angle = math.pi + args.angle[0] * pgon0.angle()/2
    else:                        # units: degrees
        turn_angle = math.radians(args.angle[0])

    try:
        (points, faces) = calc_polygons2(
            pgon0, pgon0.circumradius(), turn_angle,
            pgon1, args.ratio*pgon1.circumradius(),
            math.radians(args.angle_between_axes))
    except Exception as e:
        parser.error(e.args[0])

    if(args.x_axis_vert is not None):
        v = pgon0.N + args.x_axis_vert % pgon1.N
        P = points[v].copy()
        transl = Mat.transl(Vec(0, 0, -P[2]))
        # for i in range(len(points)):
        #    points[i][2] -= P[2]
        rot = Mat.rot_xyz(0, 0, anti_lib.angle_around_axis(
            P, Vec(1, 0, 0),  Vec(0, 0, 1)))
        points = [rot*transl*point for point in points]

    out = anti_lib.OffFile(args.outfile)
    out.print_all(points, faces)

if __name__ == "__main__":
    main()
