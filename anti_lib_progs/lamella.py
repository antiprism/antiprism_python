#!/usr/bin/env python3

# Copyright (c) 2014-2017 Adrian Rossiter <adrian@antiprism.com>
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
Create a lamella dome with a constant kite angle.
'''

import argparse
import sys
from math import cos, sin, tan, sqrt, atan, acos, pi
import anti_lib
from anti_lib import Vec

epsilon = 1e-13


def get_default_angle(pgon, radius_top, radius_bottom, height):
    a = pgon.angle()/2
    P1 = Vec(radius_top, 0, height)
    Q1 = Vec(radius_bottom * cos(a), radius_bottom * sin(a), 0.0)
    Q2 = Vec(radius_bottom * cos(-a), radius_bottom * sin(-a), 0.0)
    v0 = (Q1-P1).unit()
    v1 = (Q2-P1).unit()
    cos_a = anti_lib.safe_for_trig(Vec.dot(v0, v1))
    ang = acos(cos_a)
    return pi-ang


def get_axis_intersect(P0, P1, P2):
    mid = (P0 + P1)/2
    A = sqrt(P2[0]**2 + P2[1]**2)
    B = sqrt(mid[0]**2 + mid[1]**2)
    z_diff = mid[2] - P2[2]
    apex = P2[2] - A * z_diff / (B - A)
    return anti_lib.Vec(0, 0, apex)


def get_next_radius_height(vals0, vals1, ang, k_ang):
    R0 = vals0[0]
    h0 = vals0[1]
    R1 = vals1[0]
    h1 = vals1[1]
    h = h1 - h0
    r1 = R1*cos(ang)
    d = R0 - r1
    l = sqrt(h*h + d*d)
    e = R1*sin(ang)
    alpha = atan(e/l)           # half lower kite angle
    beta = pi - k_ang - alpha   # half higher kite angle
    h2 = h*(tan(alpha)/tan(beta)+1)
    d1 = d*(tan(alpha)/tan(beta)+1)
    R2 = R0 - d1

    return R2, h0+h2


def make_lamella_dome(pgon, radius_top, height_top, radius_bottom,
                      height_bottom, kite_ang, first_level, last_level, caps):
    N = pgon.N
    a = pgon.angle()/2
    vals = [(radius_bottom, height_bottom), (radius_top, height_top)]

    if first_level > last_level:
        first_level, last_level = last_level, first_level

    max_above = last_level
    max_below = first_level

    for n in range(max_above):
        vals.append(get_next_radius_height(vals[-2], vals[-1], a, kite_ang))
        if vals[-1][0] < 0:
            break

    vals2 = [vals[1], vals[0]]
    for n in range(-max_below):
        vals2.append(get_next_radius_height(vals2[-2], vals2[-1], a, kite_ang))
        if vals2[-1][0] < 0:
            break

    all_vals = list(reversed(vals2[2:]))+vals
    if(max_below > 0):
        all_vals = all_vals[max_below:-1]
    if(max_above < 0):
        all_vals = all_vals[0:max_above]

    points = []
    faces = []
    for n in range(len(all_vals)):
        R1 = all_vals[n][0]
        h1 = all_vals[n][1]
        not_last = (n < len(all_vals)-2)

        a_off = a * (n % 2)
        for i in range(N):
            ang = a_off + i * 2 * pi / N
            pt = Vec(R1*cos(ang), R1*sin(ang), h1)
            points.append(pt)

            if not_last:
                faces.append([n*N+i, (n+1)*N+(i+N-(n+1) % 2) % N,
                             (n+2)*N+i, (n+1)*N+(i+N+n % 2) % N])

    for cap in range(2):
        last = len(all_vals) - 1
        if caps[cap] == 'a':
            apex_index = len(points)
            if(cap == 0):
                apex = get_axis_intersect(points[0], points[1], points[N])
            else:
                apex = get_axis_intersect(points[last*N], points[last*N+1],
                                          points[(last-1)*N])
            points.append(apex)
            for i in range(N):
                if cap == 0:
                    faces.append([apex_index, i, N+i, (i+1) % N])
                else:
                    faces.append([(last-1)*N+(N + i - 1 + last % 2) % N,
                                  last*N+(N+i-1) % N,
                                  apex_index, last*N+i])
        else:
            for i in range(N):
                if cap == 0:
                    faces.append([i, N+i, (i+1) % N])
                else:
                    faces.append([(last-1)*N+(N + i - 1 + last % 2) % N,
                                  last*N+(N+i-1) % N, last*N+i])
            faces.append([(last*cap)*N + i for i in range(N)])

    return points, faces


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'polygon_fraction',
        help='number of sides of the base polygon (N), '
        'or a fraction for star polygons (N/D) (default: 12)',
        default='12',
        nargs='?',
        type=anti_lib.read_polygon)
    parser.add_argument(
        'last_level',
        help='number of levels above (+/-) to finish (default: 100,0)',
        nargs='?',
        type=int,
        default=100)
    parser.add_argument(
        'first_level',
        help='number of levels above (+/-) to start (default: 0)',
        nargs='?',
        type=int,
        default=0)
    parser.add_argument(
        '-k', '--kite-angle',
        help='side angle of kite (default: calculated)',
        type=float,
        nargs='?',
        default=-1.0)
    parser.add_argument(
        '-l', '--height',
        help='height of base antiprism (default: calculated)',
        type=float,
        nargs='?',
        default=-1.0)
    parser.add_argument(
        '-a', '--alignment',
        help='alignment: p - prism, a - antiprism (default)',
        choices=['p', 'a'],
        default='a')
    parser.add_argument(
        '-c', '--caps',
        help='cap types for bottom and top (in that order), two letters '
             'from p (polygon), and a (apex) (default: pa)',
        choices=['aa', 'ap', 'pa', 'pp'],
        default='pa')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    pgon = args.polygon_fraction
    if pgon.D > 1:
        parser.error('star polygons not currently implemented')
    pgon_a = pgon.angle() / 2

    radius_bottom = 1.0
    if args.alignment == 'a':
        radius_top = radius_bottom
        if args.height > 0:
            height_top = args.height / 2
        else:
            height_top = 0.5 * sin(pgon_a) * sqrt(2 - 1 / (cos(pgon_a / 2))**2)
        height_bottom = -height_top
    else:
        radius_top = radius_bottom * cos(pgon_a)
        if args.height > 0:
            height_top = args.height
        else:
            height_top = radius_bottom * sin(pgon_a)
        height_bottom = 0.0

    if args.kite_angle >= 0.0:
        ang = args.kite_angle*pi/180
    else:
        ang = get_default_angle(pgon, radius_top, radius_bottom,
                                height_top-height_bottom)
        print("calculated angle:"+str(180/pi*ang), file=sys.stderr)

    points, faces = make_lamella_dome(
            pgon, radius_top, height_top, radius_bottom, height_bottom,
            ang, args.first_level, args.last_level, args.caps)

    out = anti_lib.OffFile(args.outfile)
    out.print_all_pgon(points, faces, pgon)

if __name__ == "__main__":
    main()
