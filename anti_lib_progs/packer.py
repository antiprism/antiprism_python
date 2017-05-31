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
Pack balls in a sphere. The pack is seeded with two or more balls, then
subsequent balls are added one at a time in three point contact in
positions chosen by the packing method.
'''

import argparse
import math
import sys
from functools import cmp_to_key
import anti_lib
from anti_lib import Vec


# return the coordinate and radius of the circle produced by the
# intersection of two spheres with coors p1, p2 and radii r0, r1
def sphere_intersection(p0, r0, p1, r1):
    p0_p1 = p1 - p0
    len_p0_p1 = p0_p1.mag()
    # do points have the same coordinates?
    if len_p0_p1 == 0:
        return None, None

    # take a cross section through the spheres, and consider triangle
    # between sphere centres and one intersection point. compare areas by
    # herons law and base*ht/2 to find ht (ht = R = circle radius)

    sp = (r0 + r1 + len_p0_p1) / 2.0   # semi perimeter

    # do spheres intersect?
    if len_p0_p1 > sp or r0 > sp or r1 > sp:
        return None, None

    area = math.sqrt(sp*(sp - r0)*(sp - r1)*(sp - len_p0_p1))
    R = 2*area / len_p0_p1

    # find length of vec in direction p0 to p1 by cos / cos rule, hence coord
    P = p0 + p0_p1*(r0**2 + len_p0_p1**2 - r1**2)/(2*len_p0_p1**2)

    return P, R


# return coordinates of the 2 (or less) spheres of radius R that touch
# three spheres with coords p1, p2, p3 and radii r1, r2, r3
def touching_spheres(R, p0, r0, p1, r1, p2, r2):
    # Add R to the radii, the intersection points of the three sheres are
    # the centres of the touching ball(s)
    r0 += R
    r1 += R
    r2 += R

    # find centre of circle nd readius where first two spheres intersect
    pc1, rc1 = sphere_intersection(p0, r0, p1, r1)
    if pc1 is None:
        return None, None

    # the second circle is defined by th plane of first circle intersecting
    # the third sphere
    # find centre of second circle
    p2_pc1 = pc1 - p2
    p0_p1 = p1 - p0
    unit_p0_p1 = p0_p1.unit()
    len_p2_pc2 = Vec.dot(p2_pc1, unit_p0_p1)

    # is circle outside of third sphere
    if not (r2 > len_p2_pc2 > -r2):
        return None, None

    p2_pc2 = unit_p0_p1 * len_p2_pc2
    pc2 = p2 + p2_pc2

    rc2 = math.sqrt(r2**2 - len_p2_pc2**2)

    # find pt - intersection of plane of sphere centres and line joining
    # touching spheres, and R distance pt to a touching sphere centre.
    pt, R = sphere_intersection(pc1, rc1, pc2, rc2)
    if pt is None:
        return None, None

    # make a vector length R perpendicular to the plane of the sphere centres
    R_norm = Vec.cross(p1 - p0, p2 - p0)
    len_R_norm = R_norm.mag()
    if not len_R_norm:
        return None, None

    R_norm = R_norm * R/len_R_norm

    return pt + R_norm, pt - R_norm

# ------- End of Pocket functions -----------

epsilon = 1e-12


# compare function to sort on distance to to centre
def cmp_from_orig(a, b):
    a2 = a[0]**2 + a[1]**2 + a[2]**2
    b2 = b[0]**2 + b[1]**2 + b[2]**2
    if a2 < b2 - epsilon:
        return -1
    if a2 > b2 + epsilon:
        return 1
    return 0


# compare function to sort on distance to to centre
def cmp_from_orig_up(a, b):
    a2 = a[0]**2 + a[1]**2 + a[2]**2
    b2 = b[0]**2 + b[1]**2 + b[2]**2
    if a2 < b2 - epsilon:
        return -1
    if a2 > b2 + epsilon:
        return 1

    # a==b
    if a[2] < b[2] - epsilon:
        return -1
    if a[2] > b[2] + epsilon:
        return 1

    return 0


# compare function to sort on distance from container
def cmp_from_cont(a, b):
    a2 = a[0]**2 + a[1]**2 + a[2]**2
    b2 = b[0]**2 + b[1]**2 + b[2]**2
    if a2 > b2 + epsilon:
        return -1
    if a2 < b2 - epsilon:
        return 1
    return 0


# compare function to sort on distance from container then z
def cmp_from_cont_up(a, b):
    a2 = a[0]**2 + a[1]**2 + a[2]**2
    b2 = b[0]**2 + b[1]**2 + b[2]**2
    if a2 > b2 + epsilon:
        return -1
    if a2 < b2 - epsilon:
        return 1

    # a==b
    if a[2] < b[2] - epsilon:
        return -1
    if a[2] > b[2] + epsilon:
        return 1

    return 0


# compare function to sort on z
def cmp_z(a, b):
    if a[2] < b[2] - epsilon:
        return -1
    if a[2] > b[2] + epsilon:
        return 1
    return 0


# choose next ball position from pocket list
def next_ball_pos(points, pkts, r):
    for i in range(len(pkts)):
        # would a ball in this pocket overlap with any other ball
        for n in range(1, len(points)):
            overlap = 0
            v = pkts[i] - points[n]
            if v.mag2() < 4*r**2 - epsilon:
                overlap = 1
                break
        if not overlap:
            break

    if overlap:
        del pkts
        return None

    next_p = pkts[i]     # first suitable pocket
    del pkts[0:i + 1]    # remove this, and earlier unsuitable pockets
    return next_p


def find_pockets(points, new_p, r, R):
    new_pkts = []
    l = len(points)
    for i in range(l):
        # set radius of first ball (0 is the container ball)
        if i:
            ri = r
        else:
            ri = R - 2*r
        for j in range(i + 1, l):
            pt = touching_spheres(r, points[i], ri, points[j], r, new_p, r)
            if pt[0] is None:
                continue
            for p in pt:
                # check it is within container
                if p.mag2() < (R - r)**2 + epsilon:
                    new_pkts.append(p)
    return new_pkts


def positive_float(val_str):
    try:
        val = float(val_str)
    except:
        raise argparse.ArgumentTypeError('not a number')

    if val <= 0.0:
        raise argparse.ArgumentTypeError('not a positive number')

    return val


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'ball_radius',
        help='radius of balls to pack in container '
             '(default: 0.3)',
        type=anti_lib.read_positive_float,
        nargs='?',
        default=0.3)
    parser.add_argument(
        'container_radius',
        help='radius of container sphere (default: 1.0)',
        type=anti_lib.read_positive_float,
        nargs='?',
        default=1.0)
    parser.add_argument(
        '-m', '--method',
        help='packing method: up - bottom up (default), '
             'in - outside to centre, out - centre to outside, '
             'inup - outside in bottom first, '
             'outup - centre to outside bottom first.',
        choices=['up', 'in', 'out', 'inup', 'outup'],
        default='out')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    r = args.ball_radius
    R = args.container_radius
    if r > R:
        parser.error('ball_radius is greater than container_radius\n')

    points = []
    points.append(Vec(0, 0, 0))        # container sphere/ball

    if args.method == 'up':
        cmp_func = cmp_z
        points.append(Vec(0, 0, -(R - r)))  # first ball at min z
        # dummy ball point to find starting pockets
        dummy_p = Vec(r, 0, -(R - r))
    elif args.method == 'inup':
        cmp_func = cmp_from_cont_up
        points.append(Vec(0, 0, -(R - r)))  # first ball at min z
        # dummy ball point to find starting pockets
        dummy_p = Vec(r, 0, -(R - r))
    elif args.method == 'in':
        cmp_func = cmp_from_cont
        points.append(Vec(0, 0, -(R - r)))  # first ball at min z
        # dummy ball point to find starting pockets
        dummy_p = Vec(r, 0, -(R - r))
    elif args.method == 'outup':
        cmp_func = cmp_from_orig_up
        points.append(Vec(0, 0, 0))     # ball at centre
        points.append(Vec(0, 0, -2*r))    # ball below
        # dummy ball point to find starting pockets
        dummy_p = Vec(2*r, 0, 0)
    else:  # args.method=="out"
        cmp_func = cmp_from_orig
        points.append(Vec(0, 0, 0))     # ball at centre
        points.append(Vec(0, 0, 2*r))    # ball above
        # dummy ball point to find starting pockets
        dummy_p = Vec(2*r, 0, 0)

    pkts = find_pockets(points, dummy_p, r, R)

    if pkts:
        while 1:
            new_p = next_ball_pos(points, pkts, r)
            if new_p is None:
                break

            new_pkts = find_pockets(points, new_p, r, R)
            if new_pkts:
                pkts.extend(new_pkts)
                pkts.sort(key=cmp_to_key(cmp_func))
            points.append(new_p)

    print('packed {} balls of radius {} in container of radius {}\n'.format(
          len(points) - 1, r, R), file=sys.stderr)

    out = anti_lib.OffFile(args.outfile)
    out.print_all(points[1:], [])

if __name__ == "__main__":
    main()
