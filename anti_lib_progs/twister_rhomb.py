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
Twist polygons, of the same type, placed on certain fixed axes and
joined by vertices.

'''

import argparse
import sys
import os
import math
import anti_lib
from anti_lib import Vec, Mat

epsilon = 1e-13


def make_arc(v0, v1, num_segs, start_idx):
    axis = Vec.cross(v0, v1).unit()
    ang = anti_lib.angle_around_axis(v0, v1, axis)
    points = [v0]
    faces = []
    mat = Mat.rot_axis_ang(axis, ang/num_segs)
    for i in range(num_segs):
        # accumulated error
        points.append(mat * points[-1])
        faces.append([start_idx + i, start_idx + i+1])
    return points, faces


def make_frame(frame_elems, pgon, axis_angle, num_segs):
    points = []
    faces = []
    if frame_elems:
        v0 = Vec(0, 0, 1)
        v1 = Vec(-math.sin(axis_angle), 0, math.cos(axis_angle))
        v2 = v1.rot_z(pgon.angle()/2)
        v2[2] *= -1
        if 'r' in frame_elems:
            ps, fs = make_arc(v0, v1, num_segs, 0)
            points += ps
            faces += fs
            ps, fs = make_arc(v1, v2, num_segs, num_segs+1)
            points += ps
            faces += fs
        if 'a' in frame_elems:
            faces += [[len(points)+i, len(points)+i+1] for i in range(0, 6, 2)]
            points += [v0, -v0, v1, -v1, v2, -v2]

        rad = calc_polygons(pgon, 0, axis_angle, -1)[0][0].mag()
        points = [rad * p for p in points]
        faces += [[i] for i in range(len(points))]

    return points, faces


def calc_polygons(pgon, ang, angle_between_axes, sign_flag=1):
    # rotate initial vertex of first polygon by ang. Common
    # point lies on a line through vertex in the direction
    # of first axis (z-axis). Common point lies on a cylinder
    # with the circumradius for radius and second axis for axis.

    # rotate model so axis1 is on z-axis, to simplify the
    # calculation of the intersection
    rot = Mat.rot_axis_ang(Vec(0, 1, 0), angle_between_axes)

    R = pgon.circumradius()
    V = Vec(R, 0, 0).rot_z(ang)     # initial vertex turned on axis
    Q = rot * V                     # vertex in rotated model
    u = rot * Vec(0, 0, 1)          # direction of line in rotated model

    # equation of line is P = Q + tu for components x y z and parameter t
    # equation of cylinder at z=0 is x^2 + y^2 = r1^2
    a = u[0]**2 + u[1]**2
    b = 2*(Q[0]*u[0] + Q[1]*u[1])
    c = Q[0]**2 + Q[1]**2 - R**2

    disc = b**2 - 4*a*c
    if disc < -epsilon:
        raise Exception("model is not geometrically constructible")
    elif disc < 0:
        disc = 0

    # The sign flag, which changes for the range 90 to 270 degrees, allows
    # the model to reverse, otherwise the model breaks apart in this range.
    t = (-b + sign_flag*math.sqrt(disc))/(2*a)

    P = V + Vec(0, 0, t)   # the common point

    points = pgon.get_points(P)
    faces = pgon.get_faces()

    Q = rot * P
    rot_inv = Mat.rot_axis_ang(Vec(0, 1, 0), -angle_between_axes)
    points += [rot_inv*p for p in pgon.get_points(Q)]
    faces += pgon.get_faces(pgon.N*pgon.parts)

    return points, faces


def rot_reflect_pair(points, pgon, d, rev=False):
    mat0 = Mat.rot_axis_ang(Vec(0, 0, 1), (d+0.5)*pgon.angle())
    pts = [mat0 * Vec(p[0], p[1], -p[2]) for p in points]
    mat1 = Mat.rot_axis_ang(Vec(0, 0, 1), d*pgon.angle())
    if rev:
        return [mat1 * p for p in points]+pts
    else:
        return pts + [mat1 * p for p in points]


def frame_type(arg):
    if arg.strip('ra'):
        raise argparse.ArgumentTypeError(
            'frame type contains letters other than r, a')
    return arg


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
        'polygon',
        help='number of sides of polygon (default: 7) '
        '(or may be a polygon fraction, e.g. 5/2)',
        type=anti_lib.read_polygon,
        nargs='?',
        default='7')
    parser.add_argument(
        'turn_angle',
        help='amount to turn polygon on axis0 in degrees '
             '(default (0.0), or a value followed by \'e\', '
             'where 1.0e is half the central angle of an edge, '
             'which produces an edge-connected model (negative '
             'values may have to be specified as, e.g. -a=-1.0e),'
             'or a value followed by x, which is like e but with '
             'a half turn offset',
        type=read_turn_angle,
        nargs='?',
        default='0')
    parser.add_argument(
        '-n', '--number-faces',
        help='number of faces in output (default: all): '
             '0 - none (frame only), '
             '2 - cap and adjoining polygon'
             '4 - two caps and two adjoining connected polygons',
        type=int,
        choices=[0, 2, 4],
        default=-1)
    parser.add_argument(
        '-F', '--frame',
        help='include frame elements in output, any from: '
             'r - rhombic tiling edges, '
             'a - rotation axes (default: no elements)',
        type=frame_type,
        default='')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    pgon = args.polygon
    N = pgon.N
    D = pgon.D
    parts = pgon.parts
    # if N % 2 == 0:
    #    parser.error('polygon: %sfraction numerator must be odd' %
    #                 ('reduced ' if parts > 1 else ''))
    if D % 2 == 0:
        print(os.path.basename(__file__)+': warning: '
              'polygon: %sfraction denominator should be odd, '
              'model will only connect correctly at certain twist angles' %
              ('reduced ' if parts > 1 else ''),
              file=sys.stderr)

    if abs(N/D) < 3/2:
        parser.error('polygon: the polygon fraction cannot be less than 3/2 '
                     '(base rhombic tiling is not constructible)')

    axis_angle = math.acos(1/math.tan(
        math.pi*D/N)/math.tan(math.pi*(N-D)/(2*N)))
    if args.turn_angle[1] == 'e':     # units: half edge central angle
        turn_angle = args.turn_angle[0] * pgon.angle()/2
    elif args.turn_angle[1] == 'x':   # units: half edge central angle
        turn_angle = math.pi + args.turn_angle[0] * pgon.angle()/2
    else:                        # units: degrees
        turn_angle = math.radians(args.turn_angle[0])

    turn_angle_test_val = abs(math.fmod(abs(turn_angle), 2*math.pi) - math.pi)
    sign_flag = 1 - 2*(turn_angle_test_val > math.pi/2)
    try:
        (points, faces) = calc_polygons(pgon, turn_angle, axis_angle,
                                        sign_flag)
    except Exception as e:
        parser.error(e.args[0])

    if args.number_faces < 0:
        num_twist_faces = (2*N + 2)*parts
    else:
        num_twist_faces = args.number_faces*parts
    num_twist_points = num_twist_faces * N

    frame_points, frame_faces = make_frame(args.frame, pgon, axis_angle, 10)

    if num_twist_points + len(frame_points) == 0:
        parser.error('no output specified, use -f with -n 0 to output a frame')

    if D % 2 == 0:
        mat = Mat.rot_axis_ang(Vec(0, 0, 1), math.pi/N)
        points = [mat * p for p in points]

    out = anti_lib.OffFile(args.outfile)
    out.print_header(num_twist_points+2*N*len(frame_points),
                     num_twist_faces+2*N*len(frame_faces))
    if args.number_faces == -1:
        out.print_verts(rot_reflect_pair(points[:N*parts], pgon, 0, True))
        for i in range(N):
            out.print_verts(rot_reflect_pair(points[N*parts:], pgon, i))
    elif args.number_faces == 2:
        out.print_verts(points)
    elif args.number_faces == 4:
        out.print_verts(rot_reflect_pair(points[:N*parts], pgon, 0, True),)
        out.print_verts(rot_reflect_pair(points[N*parts:], pgon, 0))

    for i in range(N):
        out.print_verts(rot_reflect_pair(frame_points, pgon, i))

    for i in range(num_twist_faces):
        out.print_face(faces[0], i*N, (i//parts) % 2)

    for i in range(N):
        cur_num_points = num_twist_points + 2*i*len(frame_points)
        for j in [0, len(frame_points)]:
            out.print_faces(frame_faces, cur_num_points+j, col=2)

if __name__ == "__main__":
    main()
