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
Make a Jacoby-style dome, as described in
http://www.google.com/patents/US7900405 . Project a tiling
of unit-edged triangles, squares or crossed squares (unit edges),
at a specified height, onto a unit hemisphere, by gnomonic,
stereographic or general point projection
'''

import argparse
import sys
import math
import anti_lib
from anti_lib import Vec, Mat

epsilon = 1e-13


def make_sq_tiling(freq):
    points = []
    faces = []
    grads = 2*freq+1
    for i in range(grads):
        for j in range(grads):
            points.append(Vec(i-freq, j-freq, 0))
            if i < grads-1 and j < grads-1:
                faces.append([j + i*grads, j+1 + i*grads,
                              j+1 + (i+1)*grads, j + (i+1)*grads])
    return points, faces, points[0].mag()


def make_sq_x_tiling(freq):
    points = []
    faces = []
    grads = 2*freq+1
    for i in range(grads):
        for j in range(grads):
            points.append(Vec(i-freq, j-freq, 0))
    for i in range(grads-1):
        for j in range(grads-1):
            cent_idx = len(points)
            points.append(Vec(i-freq+0.5, j-freq+0.5, 0))
            face = [j + i*grads, j+1 + i*grads,
                    j+1 + (i+1)*grads, j + (i+1)*grads]
            for v in range(4):
                faces.append([face[v], face[(v+1) % 4], cent_idx])

    return points, faces, points[0].mag()


def get_tri_faces(t, freq):
    line0 = [t*freq+i+1 for i in range(freq)]
    line_adj = [((t+1) % 6)*freq+i+1 for i in range(freq)]
    inner_cnt = freq*(freq-1)//2
    lines = [[0]]
    n = 1 + freq*6 + t*inner_cnt
    for i in range(0, freq):
        line = [line0[i]]
        for j in range(0, i):
            line.append(n)
            n = n+1
        line.append(line_adj[i])
        lines.append(line)
    faces = []
    for i in range(freq):
        for j in range(i+1):
            if i != freq-1:
                faces.append([lines[i+1][j+1], lines[i+1][j], lines[i+2][j+1]])
            faces.append([lines[i+1][j], lines[i+1][j+1], lines[i][j]])

    return faces


def make_hexagonal_tiling(freq):
    grads = freq + 1
    points = [[0, 0, 0]]
    faces = []
    for i in range(6):
        rot = Mat.rot_axis_ang(Vec(0, 0, 1), i*math.pi/3)
        for j in range(1, grads):
            points.append(rot * Vec(j, 0, 0))

    for t in range(6):
        p1 = ((t+1) % 6)*freq + 1
        p2 = t*freq + 1
        vec = points[p1] - points[p2]
        for i in range(1, freq):
            for j in range(0, i):
                points.append(points[t*freq+1+i] + vec*(j+1))
        faces += get_tri_faces(t, freq)

    return points, faces, points[freq].mag()


def project_onto_sphere(points, ht, proj_ht):
    new_points = []
    for point in points:
        a = point[0]
        b = point[1]
        d = ht-proj_ht
        if abs(d) < epsilon:
            z = 1
            x = 0
            y = 0
        else:
            z = (a**2 + b**2 - d**2) / (a**2 + b**2 + d**2)
            x = a * (z-1) / d
            y = b * (z-1) / d

        new_points.append(Vec(x, y, -z))
    return new_points


def parallel_project(points, R, rad):
    new_points = []
    for point in points:
        x = point[0]*rad/R
        y = point[1]*rad/R
        try:
            z = math.sqrt(1 - (x**2 + y**2))
        except:
            z = 0
        new_points.append(Vec(x, y, -z))
    return new_points


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'frequency',
        help='frequency (default: 5)',
        type=anti_lib.read_positive_int,
        nargs='?',
        default=5)
    parser.add_argument(
        'tiling_ht',
        help='tiling_height (default: 4)',
        type=float,
        nargs='?',
        default=4)
    parser.add_argument(
        '-p', '--projection-height',
        help='height of projection centre below sphere centre, '
             '0.0 - gnomonic projection (default), '
             '1.0 - stereographic projection, values with a '
             'magnitude greater than 1 may not be geometrically '
             'realisable if the tiling is wider than the sphere',
        type=float,
        default=0.0)
    parser.add_argument(
        '-t', '--tiling-type',
        help='type of tiling: t - triangles, s - squares, '
             'x - squares triangulated to centres',
        choices=['t', 's', 'x'],
        default='t')
    parser.add_argument(
        '-P', '--parallel-project-radius',
        help='set radius of tiling (<=1) and parallel project onto unit '
             'sphere',
        type=float)
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    if args.tiling_ht < epsilon:
        parser.error('tiling height cannot be very small, or negative ')
    args.tiling_ht = args.tiling_ht+1

    if args.tiling_type == 't':
        points, faces, R = make_hexagonal_tiling(args.frequency)
    elif args.tiling_type == 's':
        points, faces, R = make_sq_tiling(args.frequency)
    elif args.tiling_type == 'x':
        points, faces, R = make_sq_x_tiling(args.frequency)

    rad = args.parallel_project_radius
    print(rad, file=sys.stderr)
    if(rad is None):
        points = project_onto_sphere(
            points, args.tiling_ht, args.projection_height)
    elif(rad >= 0.0 and rad <= 1.0):
        points = parallel_project(points, R, rad)
    else:
        parser.error('option -P: radius must be between 0 and 1.0')

    out = anti_lib.OffFile(args.outfile)
    out.print_all(points, faces)

if __name__ == "__main__":
    main()
