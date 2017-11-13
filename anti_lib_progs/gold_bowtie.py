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
Create a polyhedron with axial symmetry involving a golden trapezium bow-tie
'''

import argparse
import sys
import math
from math import cos, sin, tan, sqrt
import anti_lib
from anti_lib import Vec

epsilon = 1e-13
phi = (sqrt(5) + 1) / 2
gold_trap_ht = sqrt(5*phi + 2) / 2
gold_trap_diag = sqrt(2*phi + 1)


def get_belt_data(pgon, a, b):
    """Get belt construction values"""
    ang = pgon.angle()/2
    tan_a = a*sin(ang)/(a*cos(ang)+b)     # tan of half of edge a
    r_a = a/(2*tan_a)                     # radius to middle of edge a
    R_a = math.sqrt((a/2)**2 + r_a**2)    # circumradius
    return tan_a, r_a, R_a


def cup_model(pgon, model_type):
    """Construct cupula model with golden trapziums"""
    ang = pgon.angle()/2
    if model_type in ['0']:
        len_a, len_b, len_c = (phi-1)/sin(math.pi/2 - ang), 1, phi
    elif model_type in ['1', '2']:
        len_a, len_b, len_c = phi, 1, phi
    elif model_type in ['3']:
        len_a, len_b, len_c = phi, phi, 1
    ang = pgon.angle()/2
    tan_a, r_a, R_a = get_belt_data(pgon, len_a, len_b)
    P = Vec(r_a, -len_a/2, 0)
    Q = P + Vec(0, len_a, 0)
    points = []
    for i in range(pgon.N):
        points.append(P.rot_z(2*i*ang))
        points.append(Q.rot_z(2*i*ang))

    tan_b, r_b, R_b = get_belt_data(pgon, len_b, len_a)
    r_c = len_c/(2*tan(ang))
    ht2 = gold_trap_ht**2 - (r_b - r_c)**2
    #  print(ht2, file=sys.stderr)
    if ht2 < -epsilon:
        raise ValueError('model is not constructible (height)')
    elif ht2 < 0:
        ht2 = 0
    ht = math.sqrt(ht2)
    R_c = len_c/(2*sin(ang))
    V = Vec(R_c, 0, ht)
    N = pgon.N
    faces = []
    if model_type not in ['3']:
        faces.append([2*N + i for i in range(N)])
    for i in range(pgon.N):
        points.append(V.rot_z(2*i*ang))
    if model_type in ['0', '1']:
        faces.append([3*N + i for i in range(N)])
        U = Vec(R_c, 0, -ht)
        for i in range(pgon.N):
            points.append(U.rot_z(2*i*ang))
    elif model_type in ['3']:
        faces.append([i for i in range(2*N)])

    for i in range(N):
        if model_type in ['0']:
            faces.append([2*i, 3*N+i, 2*i+1, 2*N+i])
        else:
            faces.append([2*i, 2*i+1, 2*N+i])
        faces.append([2*i+1, (2*i+2) % (2*N),
                      2*N + ((i+1) % N), 2*N + i])

    if model_type in ['0', '1']:
        for i in range(N):
            if model_type in ['1']:
                faces.append([2*i, 2*i+1, 3*N+i])
            faces.append([2*i+1, (2*i+2) % (2*N),
                          3*N + ((i+1) % N), 3*N + i])

    elif model_type in ['2']:
        belt_ht2 = gold_trap_ht**2 - (r_b - r_a)**2
        if belt_ht2 < -epsilon:
            raise ValueError('model is not constructible (belt height)')
        elif belt_ht2 < 0:
            belt_ht2 = 0
        belt_ht = math.sqrt(belt_ht2)
        for i, p in enumerate(points):
            points[i][2] += belt_ht/2
        points += [Vec(p[0], p[1], -p[2]).rot_z(ang) for p in points]
        faces += [[idx+3*N for idx in f] for f in faces]
        for i in range(2*N):
            faces.append([i, (i+1) % (2*N), 3*N + i, 3*N + (2*N+i-1) % (2*N)])

    elif model_type in ['3']:
        for i, p in enumerate(points):
            points[i][2] -= ht
        points += [Vec(p[0], p[1], -p[2]) for p in points[:2*N]]

        def new_v(v): v+3*N if v < 2*N else v

        for i in range(len(faces)):
            faces.append([new_v(v) for v in faces[i]])

    return points, faces


def bifrustum_model(pgon, arg_type):
    """Construct bifrustum model with golden trapziums"""
    N = pgon.N
    ang = pgon.angle()/2
    R0 = phi/(2*sin(ang))
    R1 = 1/(2*sin(ang))
    ht = sqrt(phi**2 - (R0-R1)**2)
    points = []
    points += [Vec(R0, 0, ht).rot_z(2*i*ang) for i in range(N)]
    points += [Vec(R1, 0, 0).rot_z(2*i*ang) for i in range(N)]
    points += [Vec(R0, 0, -ht).rot_z(2*i*ang) for i in range(N)]

    faces = []
    faces.append([i for i in range(N)])
    faces.append([i+2*N for i in range(N)])
    for i in range(N):
        faces.append([i, (i+1) % N, N + (i+1) % N, N + i])
        faces.append([N + i, N + (i+1) % N, 2*N + (i+1) % N, 2*N + i])
    return points, faces


def tri_antiprism_pt(pgon, a, b, c):
    """Return construct point for antiprism model"""
    s = (a + b + c) / 2
    alt = 2*math.sqrt(s*(s-a)*(s-b)*(s-c)) / a  # planar height of A
    ang = pgon.angle()/2
    R = pgon.circumradius(a)                  # polygon circumradius (side a)
    r = pgon.inradius(a)                      # polygon inradius (side a)
    a0 = math.sqrt(b**2 - alt**2)             # from C to perpendicular from A
    y = a0 - a/2
    x = math.sqrt(R**2 - y**2)
    alt_proj_x = x - r                        # projection of alt onto x dir
    z = math.sqrt(alt**2 - alt_proj_x**2) / 2
    ang_to_A = math.atan2(y, x)               # adjust point for dih axis
    ang_off = (ang_to_A - ang)/2
    return Vec(x, y, z).rot_z(-ang_off)


def antiprism_model(pgon, arg_type):
    """Construct antiprism model with golden trapziums"""
    N = pgon.N
    ang = pgon.angle()/2
    P = tri_antiprism_pt(pgon, gold_trap_diag, 1, phi)
    Q = Vec(P[0], -P[1], -P[2])
    points = []
    for i in range(N):
        points.append(P.rot_z(2*i*ang))
    for i in range(N):
        points.append(Q.rot_z(2*i*ang))
    R = points[N-1] + (P - Q)*phi
    for i in range(N):
        points.append(R.rot_z(2*i*ang))
    S = Vec(R[0], -R[1], -R[2])
    for i in range(N):
        points.append(S.rot_z(2*i*ang))

    faces = []
    faces.append([2*N + i for i in range(N)])
    faces.append([3*N + i for i in range(N)])
    for i in range(N):
        faces.append([i, 2*N + i, 2*N + ((i+1) % N)])
        faces.append([N+i, 3*N + i, 3*N + ((i-1) % N)])
        faces.append([i, 2*N + ((i+1) % N),
                      ((i+1) % N), N + ((i+1) % N)])
        faces.append([i, N + ((i+1) % N), 3*N + i, N + i])
    return points, faces


def main():
    """Entry point"""
    epilog = '''
notes:
  Depends on anti_lib.py. Not all models are constructible. Use of golden
  trapeziums in polyhedra proposed by Dave Smith -
  https://hedraweb.wordpress.com/

examples:
  Pentagonal orthobicupola type
  gold_bowtie.py 5 1 | antiview

  Petagramatic elongated gyrobicupola type
  gold_bowtie.py 5/2 2 | antiview

'''

    parser = argparse.ArgumentParser(formatter_class=anti_lib.DefFormatter,
                                     description=__doc__, epilog=epilog)

    parser.add_argument(
        'polygon_fraction',
        help='number of sides of the base polygon (N), '
             'or a fraction for star polygons (N/D) (default: 5)',
        default='5',
        nargs='?',
        type=anti_lib.read_polygon)
    parser.add_argument(
        'poly_type',
        help=''']polyhedron type (default: 0):
0 - like a square barrel polyhedron
1 - orthobicupola
2 - elongated gyrobicupola
3 - inverted bicupola
4 - inverted bifrustum
5 - like a double square barrel polyhedron''',
        choices=['0', '1', '2', '3', '4', '5'],
        default='1',
        nargs='?')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()
    pgon = args.polygon_fraction

    if args.poly_type in ['0', '1', '2', '3']:
        points, faces = cup_model(pgon, args.poly_type)
    elif args.poly_type == '4':
        points, faces = bifrustum_model(pgon, args.poly_type)
    elif args.poly_type == '5':
        points, faces = antiprism_model(pgon, args.poly_type)
    else:
        parser.error('unknown polyhedron type')

    out = anti_lib.OffFile(args.outfile)
    out.print_all_pgon(points, faces, pgon)

if __name__ == "__main__":
    main()
