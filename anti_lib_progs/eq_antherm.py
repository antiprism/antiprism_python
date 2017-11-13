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
Create an antihermaphrodite with equilateral triangles
'''

import argparse
import sys
from math import cos, sin, tan, sqrt, radians
import anti_lib
from anti_lib import Vec


def make_equ_antiherm(pgon, dih_ang):
    """Make a hermaphrodite with equilateral triangles"""
    N = pgon.N
    a = pgon.angle()/2
    R = pgon.circumradius()
    tri_alt = sqrt(3)/2
    tri_ht = tri_alt * sin(dih_ang)
    R2 = pgon.inradius() + tri_alt * cos(dih_ang)

    points = [Vec(0, 0, 0)]*(2*N+1)
    for i in range(N):
        points[i] = Vec(R*cos(2*i*a), R*sin(2*i*a), 0)
        points[i+N] = Vec(R2*cos(2*i*a+a), R2*sin(2*i*a+a), tri_ht)

    A = sqrt(points[0][0]**2 + points[0][1]**2)
    mid = anti_lib.centroid([points[N], points[2*N-1]])
    B = sqrt(mid[0]**2 + mid[1]**2)
    z_diff = mid[2] - points[0][2]
    points[2*N][2] = A * z_diff / (A - B)

    faces = []
    faces.append([i for i in range(N)])      # bottom

    for i in range(N):
        faces.append([i, (i + 1) % N, N + i])
        faces.append([N + i, 2*N, N + (i + 1) % N, (i + 1) % N])

    return points, faces


def main():
    """Entry point"""

    epilog = '''
notes:
  Depends on anti_lib.py.

examples:
  Pentagonal hermaphrodite with equilateral triangles angled at 100 degrees
  from base
  eq_antherm.py -a 100 5 | antiview
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
        '-a', '--angle',
        help='dihedral angle at base (default: 90)',
        type=float,
        default='90')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    pgon = args.polygon_fraction
    dih = radians(args.angle)

    points, faces = make_equ_antiherm(pgon, dih)

    out = anti_lib.OffFile(args.outfile)
    out.print_all_pgon(points, faces, pgon)

    c = sqrt(3)/6
    o = pgon.inradius() + c*cos(dih)
    h = o/tan(dih) + c*sin(dih)
    print("0,0,%f" % (h), file=sys.stderr)

if __name__ == "__main__":
    main()
