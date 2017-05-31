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
Create a jitterbug model for a general polygon base. The transformation
includes, if constructible, models corresponding to the antiprism,
snub-antiprisms and gyrobicupola for that base. The transformation is
controlled by the angle of an equatorial edge.
'''

import argparse
import sys
import math
from math import cos, sin, tan, sqrt
import anti_lib
from anti_lib import Vec

EPSILON = 1e-15   # limit for precision considerations


def print_antiprism_angle(p_ang):
    """Print the angle corresponding to an antiprism"""
    edge = 1                                    # polygon edge length
    rad = edge/(2*sin(p_ang))                   # circumradius of base polygon
    dist = 2*rad*sin(p_ang/2)
    msg = 'angle parameter(s) for antiprism: '
    val = 1 - dist**2
    if val < -EPSILON:
        msg += 'antiprism is not constructible'
    else:
        if p_ang > math.pi/2 + EPSILON:
            msg += '-x '
        if val < 0.0:
            val = 0.0
        ht = sqrt(val)
        side = 2*rad*sin(p_ang/2)
        msg += str(math.degrees(math.atan2(ht, side)))

    print(msg, file=sys.stderr)


def get_principal_verts(p_ang, ang, cross_flag):
    """Calculate example vertex coordinates"""
    edge = 1                                    # polygon edge length
    rad = edge/(2*sin(p_ang))                   # circumradius of base polygon
    belt_irad = edge*cos(ang)/(2*tan(p_ang/2))  # inradius of belt polygon
    tri_ht = edge * sqrt(3)/2

    # Cut polygon cylinder (radius rad) by plane normal to a belt edge
    # and the through its centre point
    # Distance to this ellipse should be tri_ht
    S2 = sin(ang)**2
    if abs(S2) < EPSILON:
        x = (1-2*cross_flag)*rad

        val = tri_ht**2 - (belt_irad - x)**2
        if val < -EPSILON:
            raise ValueError('triangle cannot reach base polygon '
                             '(2nd coordinate)')
        elif val < 0.0:
            val = 0.0

        y = sqrt(val)

    else:
        a = (1-S2)
        b = 2*belt_irad*S2
        c = S2*tri_ht**2 - rad**2 - S2*belt_irad**2

        disc = b**2 - 4*a*c
        if disc < - EPSILON:
            raise ValueError('triangle cannot reach base polygon '
                             '(1st coordinate)')
        elif disc < 0.0:
            disc = 0.0

        if not a:
            raise ValueError('triangle cannot reach base polygon '
                             '(vertical edges)')

        x = (-b + (1-2*cross_flag)*sqrt(disc))/(2*a)

        val = (rad**2 - x**2)/S2
        if val < -EPSILON:
            raise ValueError('triangle cannot reach base polygon '
                             '(2nd coordinate)')
        elif val < 0.0:
            val = 0.0

        y = sqrt(val)

    return (Vec(x, y*sin(ang), y*cos(ang)),                          # base
            Vec(belt_irad, 0.5*edge*cos(ang), 0.5*edge*sin(ang)))    # belt


def make_jitterbug(pgon, ang, cross, right_fill, left_fill, delete_main_faces):
    """Make the jitterbug model"""
    p_ang = pgon.angle()/2
    A, B = get_principal_verts(p_ang, ang, cross)

    A2 = Vec(A[0], -A[1], -A[2]).rot_z(-p_ang)
    A = Vec(A2[0], A2[1], -A2[2]).rot_z(p_ang)
    B2 = Vec(B[0], -B[1], -B[2])

    N = pgon.N
    points = []
    for pt in [A, A2, B2, B]:
        points += [pt.rot_z(i*2*p_ang) for i in range(N)]

    faces = []
    if not delete_main_faces:
        faces.append([i for i in range(N)])          # top
        faces.append([2*N-1 - i for i in range(N)])  # bottom

        # side triangles
        for i in range(N):
            faces.append([2*N+i, 3*N+i, (i+N) % N])
            faces.append([3*N+i, 2*N+(i+1) % N, N + (i+1+N) % N])

    # fill holes
    # Top hole: i, (i+1) % N, 2*N+(i+1) % N, 3*N+i
    # Bottom hole: N+i, N+(i-1+N) % N, 2*N+(i-1+N) % N, 3*N+(i-1) % N
    if right_fill:
        for i in range(N):
            faces.append([(i+1) % N, i, 3*N+i])
            faces.append([3*N+i, (i+1) % N, 2*N+(i+1) % N])
            faces.append([N+(i-1+N) % N, N+i, 2*N+(i-1+N) % N])
            faces.append([2*N+(i-1+N) % N, N+i, 3*N+(i-1) % N])
    if left_fill:
        for i in range(N):
            faces.append([i, (i+1) % N, 2*N+(i+1) % N])
            faces.append([2*N+(i+1) % N, 3*N+i, i])
            faces.append([N+i, N+(i-1+N) % N, 3*N+(i-1) % N])
            faces.append([3*N+(i-1) % N, 2*N+(i-1+N) % N, N+(i-1+N) % N])

    return points, faces


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'polygon_fraction',
        help='number of sides of the base polygon (N), '
        'or a fraction for star polygons (N/D)',
        default='3',
        type=anti_lib.read_polygon)
    parser.add_argument(
        'angle',
        help='angle to rotate a belt edge from horizontal',
        type=float,
        default=0.0)
    parser.add_argument(
        '-x', '--cross',
        help='alternative model type that will make a '
        'cross-antiprism for fractions greater than 1/2',
        action='store_true')
    parser.add_argument(
        '-l', '--left-fill',
        help='fill holes with triangles using the "right" '
        'diagonal',
        action='store_true')
    parser.add_argument(
        '-r', '--right-fill',
        help='fill holes with triangles using the "left" '
        'diagonal',
        action='store_true')
    parser.add_argument(
        '-d', '--delete-main-faces',
        help='delete base polygons and side triangles, '
        'keep all vertices and any faces specified with -l, -r',
        action='store_true')
    parser.add_argument(
        '-p', '--print-info',
        help='print information about the model to the screen '
             '(currently only parameters to make antiprism)',
        action='store_true')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    pgon = args.polygon_fraction

    if args.print_info:
        print_antiprism_angle(pgon.angle()/2)

    try:
        points, faces = make_jitterbug(
            pgon, math.radians(args.angle), args.cross, args.right_fill,
            args.left_fill, args.delete_main_faces)
    except ValueError as e:
        parser.error(e.args[0])

    out = anti_lib.OffFile(args.outfile)
    out.print_all_pgon(points, faces, pgon)

if __name__ == "__main__":
    main()
