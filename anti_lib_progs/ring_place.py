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
Place maximum radius rings of contacting balls around points on a sphere.
Input is a list of coordinates, one set per line.
'''

import argparse
import sys
import math
import re
import anti_lib
from anti_lib import Vec, Mat


def ring_ball_ang(N, a):
    target = math.sin(math.pi/(2*N))
    ang_range = [0, a/2]
    for i in range(100):
        diff = (ang_range[1] - ang_range[0])/2
        b = (ang_range[0] + ang_range[1])/2
        if math.sin(b/2)/math.sin(a-b) > target:
            ang_range[1] -= diff
        else:
            ang_range[0] += diff

    return b


def make_ring(R, N, a):
    b = ring_ball_ang(N, a)
    P = Vec(R*math.sin(a-b), 0, R*math.cos(a-b))
    return [P.rot_z(2*math.pi*i/N) for i in range(N)], b


def read_coords():
    points = []
    while 1:
        line = sys.stdin.readline()
        if line == '\n':
                continue

        if line == "":
            break

        m = re.search('^ *([^ ,]+) *,? *([^ ,]+) *,? *([^ ,\n]+) *$', line)
        if not m:
            sys.stderr.write(
                'error: did not find x, y and z values in following '
                'line (1):\n')
            sys.stderr.write(line)
            sys.exit(1)
        else:
            try:
                points.append(Vec(*[float(m.group(i)) for i in range(1, 3+1)]))
            except:
                sys.stderr.write(
                    'error: did not find x, y and z values in following '
                    'linei (2):\n')
                sys.stderr.write(line)
                sys.exit(1)

    return points


def find_minimum_separation(points):
    min_dist2 = 1e100
    for i in range(len(points)-1):
        for j in range(i+1, len(points)):
            v = points[i] - points[j]
            dist2 = v.mag2()
            if dist2 < min_dist2:
                min_dist2 = dist2

    return math.sqrt(min_dist2)


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'num_balls_on_ring',
        help='Number of balls on each ring',
        type=int,
        nargs='?',
        default=10)
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    ring_centres = read_coords()
    if not len(ring_centres):
        parser.error('no coordinates in input')

    R = ring_centres[0].mag()
    dist = find_minimum_separation(ring_centres)
    a = math.asin(dist/(2*R))
    ball_points, ball_ang = make_ring(R, args.num_balls_on_ring, a)
    print('ball radius = %.14f' % (2*R*math.sin(ball_ang/2)), file=sys.stderr)

    out = anti_lib.OffFile(args.outfile)
    out.print_header(len(ring_centres)*len(ball_points), 0)
    for cent in ring_centres:
        mat = Mat.rot_from_to(Vec(0, 0, 1), cent)
        out.print_verts([mat * p for p in ball_points])

if __name__ == "__main__":
    main()
