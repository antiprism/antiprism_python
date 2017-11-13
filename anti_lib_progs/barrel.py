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
Create a cyclic polyhedron with one or two bands
of equatorial squares, oriented like diamonds.
'''
import argparse
import sys
from math import cos, sin, tan, sqrt
import anti_lib


def make_1_band_model(pgon, trapezo):
    """Make square barrel model with one band of squares"""
    N = pgon.N
    a = pgon.angle() / 2
    R = 1 / sqrt(2) / sin(a)
    ht = sqrt(2)
    R_ht = 1 / sqrt(2) / tan(a)
    points = [anti_lib.Vec(0, 0, 0)] * 3 * N
    for i in range(N):
        points[i] = anti_lib.Vec(
            R_ht * cos(2 * i * a), R_ht * sin(2 * i * a), ht / 2)
        points[i + N] = anti_lib.Vec(R * cos(2 * i * a + a),
                                     R * sin(2 * i * a + a), 0)
        points[i + 2 * N] = anti_lib.Vec(R_ht * cos(2 * i * a),
                                         R_ht * sin(2 * i * a), -ht / 2)

    if trapezo:
        A = sqrt(points[N][0]**2 + points[N][1]**2)
        mid = anti_lib.centroid([points[0], points[1]])
        B = sqrt(mid[0]**2 + mid[1]**2)
        z_diff = mid[2] - points[N][2]
        apex = points[N][2] - A * z_diff / (B - A)
        points.append(anti_lib.Vec(0, 0, apex))
        points.append(anti_lib.Vec(0, 0, -apex))

    faces = []
    if not trapezo:
        faces.append([i for i in range(N)])      # top
        faces.append([i + 2 * N for i in range(N)])  # bottom

    for i in range(N):
        faces.append([i, i + N, i + 2 * N, (i - 1 + N) % N + N])
        if trapezo:
            faces.append([i, 3 * N, (i + 1) % N, i + N])
            faces.append([i + N, (i + 1) % N + 2 * N, 3 * N + 1, i + 2 * N])
        else:
            faces.append([i, (i + 1) % N, i + N])
            faces.append([i + N, (i + 1) % N + 2 * N, i + 2 * N])

    return points, faces


def make_2_band_model(pgon, trapezo):
    """Make square barrel model with two bands of squares"""
    N = pgon.N
    a = pgon.angle() / 2
    R = 1 / (sqrt(2) * sin(a))
    ht = sqrt(1 - 0.5 / (cos(a / 2))**2)
    points = [anti_lib.Vec(0, 0, 0)] * 4 * N
    for i in range(N):
        points[i] = anti_lib.Vec(
            R * cos(2 * i * a + a / 2), R * sin(2 * i * a + a / 2), ht / 2)
        points[i + N] = anti_lib.Vec(R * cos(2 * i * a - a / 2),
                                     R * sin(2 * i * a - a / 2), -ht / 2)

    for i in range(N):
        # Add in the fourth point of the squares
        points[i + 2 * N] = points[i] - \
            points[(i + 1) % N + N] + points[(i + 1) % N]
        points[i + 3 * N] = points[i + N] - points[i] + points[(i + 1) % N + N]

    if trapezo:
        A = sqrt(points[1][0]**2 + points[1][1]**2)
        mid = (points[2 * N] + points[2 * N + 1]) / 2
        B = (1 - 2 * (N < 4 * pgon.D)) * sqrt(mid[0]**2 + mid[1]**2)
        z_diff = mid[2] - points[1][2]
        apex = points[1][2] - A * z_diff / (B - A)
        points.append(anti_lib.Vec(0, 0, apex))
        points.append(anti_lib.Vec(0, 0, -apex))

    faces = []
    if not trapezo:
        faces.append([i + 2 * N for i in range(N)])  # top
        faces.append([i + 3 * N for i in range(N)])  # bottom

    for i in range(N):
        faces.append([i, (i + 1) % N + N, (i + 1) % N, i + 2 * N])
        faces.append([i + N, i, (i + 1) % N + N, i + 3 * N])
        if trapezo:
            faces.append([i + 3 * N, (i + 1) % N + N, (i + 1) %
                          N + 3 * N, 4 * N + 1])
            faces.append([i + 2 * N, (i + 1) % N, (i + 1) % N + 2 * N, 4 * N])
        else:
            faces.append([i + 3 * N, (i + 1) % N + N, (i + 1) % N + 3 * N])
            faces.append([i + 2 * N, (i + 1) % N, (i + 1) % N + 2 * N])

    return points, faces


def main():
    """Entry point"""
    epilog = '''
notes:
  Depends on anti_lib.py.

examples:
  Cuboctahedron
  barrel.py 4 | antiview

  Pentagonal square barrel capped with pyramids ("trapezobarrel")
  barrel.py -t 5 | antiview

  Pentagonal square barrel with two bands of squares
  barrel.py -n 2 5 | antiview
'''

    parser = argparse.ArgumentParser(formatter_class=anti_lib.DefFormatter,
                                     description=__doc__, epilog=epilog)

    parser.add_argument(
        'polygon_fraction',
        help='number of sides of the base polygon (N), '
        'or a fraction for star polygons (N/D) (default: 6)',
        default='6',
        nargs='?',
        type=anti_lib.read_polygon)
    parser.add_argument(
        '-n', '--number-bands',
        help='number of bands of squares (default: 1)',
        choices=['1', '2'],
        default='1')
    parser.add_argument(
        '-t', '--trapezo',
        help='make a "trapezobarrel" by extending the triangles '
             'into kites (like a trapezohedron can be made from '
             'an antiprism)',
        action='store_true')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    pgon = args.polygon_fraction
    if args.number_bands == '1':
        points, faces = make_1_band_model(pgon, args.trapezo)
    else:
        if pgon.N < 2 * pgon.D:
            parser.error('polyhedron is not constructible (fraction > 1/2)')
        points, faces = make_2_band_model(pgon, args.trapezo)

    out = anti_lib.OffFile(args.outfile)
    out.print_all_pgon(points, faces, pgon)


if __name__ == "__main__":
    main()
