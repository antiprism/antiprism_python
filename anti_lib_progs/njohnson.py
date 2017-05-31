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
Create a Johnson-based model, from J88, J89, J90, with a general polygon base.
'''

import argparse
import sys
import math
from math import cos, sin, tan, acos, atan2, sqrt
import anti_lib
from anti_lib import Vec, Mat

epsilon = 1e-15
edge = 1


def get_three_line_vert(A0, A1, B0, B1, ridge_down):
    A_edge = A1 - A0
    B_edge = B0 - A0
    cos_gamma = Vec.dot(A_edge.unit(), B_edge.unit())
    if abs(cos_gamma) > 1:
        cos_gamma /= abs(cos_gamma)
    gamma = acos(cos_gamma)
    d = (edge/2)*tan(gamma/2)
    h2 = 3*edge/4 - d**2
    if h2 < -epsilon:
        raise ValueError(
            'Not possible to calculate a three-line-fill triangle')
    h = (1 - 2*ridge_down) * sqrt(h2)

    A_mid = (A0 + A1) / 2
    B_mid = (B0 + B1) / 2
    mid_dir = (B_mid - A_mid).unit()
    norm = Vec.cross(B_edge, A_edge).unit()
    P_from_A_mid = mid_dir*d + norm*h
    return A_mid + P_from_A_mid


def get_pstt_cap_verts(pgon, ang, ridge_down):
    p_ang = pgon.angle()/2
    rad = pgon.circumradius(edge)

    # Polygon vertex
    A = Vec(rad, 0, 0)

    # First and last square vertices
    B = Vec(rad + edge*cos(ang)*cos(p_ang),
            edge*cos(ang)*sin(p_ang),
            edge*sin(ang))

    B2 = Vec(B[0], -B[1], B[2]).rot_z(2*p_ang)

    A_minus1 = A.rot_z(-2*p_ang)
    B2_minus1 = B2.rot_z(-4*p_ang)
    C = get_three_line_vert(A, A_minus1, B, B2_minus1, ridge_down)

    return A, B, B2, C


def get_pstts_cap_verts(pgon, ang, ridge_down):
    p_ang = pgon.angle()/2
    rad = pgon.circumradius(edge)

    # Polygon vertex
    A = Vec(rad, 0, 0)

    # First and last square vertices
    B = Vec(rad + edge*cos(ang)*cos(p_ang),
            edge*cos(ang)*sin(p_ang),
            edge*sin(ang))

    B2 = Vec(B[0], -B[1], B[2])

    # Triangle vertex
    C = Vec(0, 0, 0)
    mid = Vec(B[0], 0, B[2])  # M
    s1 = Vec(B[0]-rad, 0, B[2]).mag()   # length from  A to mid of B, B2

    # print("C = %g, sqrt(3)/(2*sin(alpha)) = %g\n" %
    #      (cos(ang), sqrt(3/4)/sin(p_ang)), file=sys.stderr)

    # b is angle CAM
    cos_b = edge/(2*s1)  # triangle pair is part of a disphenoid
    if abs(cos_b) > 1+epsilon:
        raise ValueError('Not possible to calculate ridge vertex')
    elif abs(cos_b) >= 1:
        b = 0
    else:
        b = acos(cos_b)

    mid_ang = atan2(mid[2], mid[0]-rad)

    if ridge_down:
        C_ang = -b + mid_ang   # down
    else:
        C_ang = b + mid_ang    # up

    C = Vec(rad + edge*cos(C_ang), 0, edge*sin(C_ang))
    return A, B, B2, C


def j88_get_principal_verts(pgon, ang, flags):
    bad = flags.strip('AB')
    if bad:
        raise ValueError('Unrecognised form flag(s) \''+bad+'\'')
    ridge = 'A' not in flags
    belt_back = 'B' in flags
    A, B, B2, C = get_pstts_cap_verts(pgon, ang, ridge)

    p_ang = pgon.angle()/2
    sq_mid = (B + B2.rot_z(2*p_ang)) / 2
    sq_mid_rad = Vec(sq_mid[0], sq_mid[1], 0).mag()
    A_rad = A[0]
    tri_ht2 = 3*edge/4 - (A_rad - sq_mid_rad)**2
    if tri_ht2 < -epsilon:
        raise ValueError(
            'Not possible to calculate upper equilateral triangle')

    A2_ht = (1 - 2*belt_back)*sqrt(tri_ht2) + B[2]
    A2 = Vec(A_rad, 0, A2_ht).rot_z(p_ang)

    mid_B_B1 = Vec(B[0], 0, B[2])
    ax = C - mid_B_B1
    rot = Mat.rot_axis_ang(ax, math.pi)
    pt = A - mid_B_B1
    C2 = rot * pt + mid_B_B1

    return A, B, B2, C, A2, C2


def j88_make_poly(pgon, ang, flags):
    principal_pts = j88_get_principal_verts(pgon, ang, flags)
    N = pgon.N
    p_ang = pgon.angle()/2
    points = []
    for pt in principal_pts:
        points += [pt.rot_z(i*2*p_ang) for i in range(N)]

    faces = []
    # base polygons
    faces.append([i for i in range(N)])
    faces.append([4*N + i for i in range(N)])

    # squares
    for i in range(N):
        faces.append([i, N + i, 2*N + (i + 1) % N, (i + 1) % N])

    # side triangles
    for i in range(N):
        faces.append([i, 3*N + i, N + i])
        faces.append([i, 2*N + i, 3*N + i])

    # vertical mirror triangles
    for i in range(N):
        faces.append([4*N + i, N + i, 2*N + (i+1) % N])

    for i in range(N):
        # pentagon to form pyramid with C2
        pent = [4*N + i, N + i, 3*N + i, 2*N + i, 4*N + (i+N-1) % N]
        for p in range(5):
            faces.append([5*N+i, pent[p], pent[(p+1) % 5]])

    return points, faces


def j89_get_principal_verts(pgon, ang, flags):
    bad = flags.strip('ABC')
    if bad:
        raise ValueError('Unrecognised form flag(s) \''+bad+'\'')
    ridge_down = 'A' in flags
    ridge2_down = 'B' in flags
    belt_back = 'C' in flags
    pgon.N *= 2
    # p_ang2 = pgon2.angle()
    A, B, B2, C = get_pstt_cap_verts(pgon, ang, ridge_down)
    pgon.N //= 2

    p_ang = pgon.angle()/2
    A2_rad = pgon.circumradius(edge)  # circumradius of top polygon
    sq_mid = (B + B2) / 2
    sq_mid_rad = Vec(sq_mid[0], sq_mid[1], 0).mag()
    tri_ht2 = 3/4 - (A2_rad - sq_mid_rad)**2
    if tri_ht2 < -epsilon:
        raise ValueError(
            'Not possible to calculate upper equilateral triangle')

    A2_ht = (1 - 2*belt_back)*sqrt(tri_ht2) + B[2]
    A2 = Vec(A2_rad, 0, A2_ht).rot_z(p_ang/2)

    A2_minus1 = A2.rot_z(-2*p_ang)
    B2_minus1 = B2.rot_z(-2*p_ang)
    C2 = get_three_line_vert(A2_minus1, A2, B2_minus1, B, ridge2_down)
    return [pt.rot_z(p_ang/2) for pt in [A, B, B2, C, A2, C2]]


def j89_make_poly(pgon, ang, flags):
    N = pgon.N
    p_ang = pgon.angle()/4
    principal_pts = j89_get_principal_verts(pgon, ang, flags)
    points = [principal_pts[0].rot_z(i*2*p_ang) for i in range(2*N)]
    for pt in principal_pts[1:]:
        points += [pt.rot_z(i*4*p_ang) for i in range(N)]

    faces = []
    # base polygons
    faces.append([i for i in range(2*N)])
    faces.append([5*N + i for i in range(N)])

    # squares
    for i in range(N):
        faces.append([2*i, 2*N + i, 3*N + i, (2*i+1) % (2*N)])

    # vertical mirror triangles above squares
    for i in range(N):
        faces.append([2*N + i, 3*N + i, 5*N + i])

    # bottom vertical mirror triangles
    for i in range(N):
        faces.append([2*i, (2*i-1) % (2*N), 4*N + i])

    # bottom slant triangles
    for i in range(N):
        faces.append([(2*i-1) % (2*N), 3*N + (i-1) % N, 4*N + i])
        faces.append([2*i, 4*N + i, 2*N + i])

    # top vertical mirror triangles
    for i in range(N):
        faces.append([5*N + (i-1) % N, 5*N + i, 6*N + i])

    # top slant triangles
    for i in range(N):
        faces.append([5*N + i, 2*N + i, 6*N + i])
        faces.append([5*N + (i-1) % N, 6*N + i, 3*N + (i-1) % N])

    # mid slant triangles
    for i in range(N):
        faces.append([6*N + i, 2*N + i, 4*N + i])
        faces.append([3*N + (i-1) % N, 6*N + i, 4*N + i])

    return points, faces


def j90_get_principal_verts(pgon, ang, flags):
    bad = flags.strip('AB')
    if bad:
        raise ValueError('Unrecognised form flag(s) \''+bad+'\'')
    ridge = 'A' not in flags
    belt_back = 'B' in flags
    A, B, B2, C = get_pstts_cap_verts(pgon, ang, ridge)

    p_ang = pgon.angle()/2

    # print("sin(mid_ang=%g (%g), cos(mid_ang)=%g(%g)"
    #      % (sin(mid_ang), sqrt((1-cos(ang)**2)/(1-(sin(p_ang)*cos(ang))**2)),
    #         cos(mid_ang),
    #         sqrt((cos(p_ang)*cos(ang))**2/(1-(sin(p_ang)*cos(ang))**2))),
    #      file=sys.stderr)

    # find height of equatorial triangle on vertical mirror plane
    sq_mid = (B + B2.rot_z(2*p_ang)) / 2
    # z_diff_b = vec_len([C[0], C[1], 0]) - vec_len([sq_mid[0], sq_mid[1], 0])
    z_diff = (C - sq_mid.rot_z(-p_ang))[0]
    # print("z_diff:\n%0.17g\n%0.17g\n" % (z_diff, z_diff_b), file=sys.stderr)
    tri_ht2 = edge**2*3/4 - z_diff**2
    if tri_ht2 < -epsilon:
        raise ValueError('Not possible to calculate equatorial triangle')
    elif tri_ht2 <= 0:
        tri_ht = 0
    else:
        tri_ht = sqrt(tri_ht2)

    if(belt_back):
        tri_ht *= -1

    half_ht = (B[2] + C[2] + tri_ht)/2
    A[2] -= half_ht
    B[2] -= half_ht
    B2[2] -= half_ht
    C[2] -= half_ht

    return A, B, B2, C


def j90_make_poly(pgon, ang, flags):
    principal_pts = j90_get_principal_verts(pgon, ang, flags)
    N = pgon.N
    p_ang = pgon.angle()/2
    points = []
    for s in [1, -1]:
        for pt in principal_pts:
            pt[2] *= s
            points += [pt.rot_z(i*2*p_ang + s*p_ang/2) for i in range(N)]

    faces = []
    # base polygons
    faces.append([i for i in range(N)])
    faces.append([4*N + i for i in range(N)])

    # squares
    for i in range(N):
        faces.append([i, N + i, 2*N + (i + 1) % N, (i + 1) % N])
        faces.append([4*N + i, 4*N + N + i, 4*N + 2*N + (i + 1) % N,
                      4*N + (i + 1) % N])

    # side triangles
    for i in range(N):
        faces.append([i, 3*N + i, N + i])
        faces.append([4*N + i, 4*N + 3*N + i, 4*N + N + i])
        faces.append([i, 2*N + i, 3*N + i])
        faces.append([4*N + i, 4*N + 2*N + i, 4*N + 3*N + i])

    # mirror triangles
    for i in range(N):
        faces.append([3*N + i, 5*N + i, 6*N + (i + 1) % N])
        faces.append([7*N + i, N + (i-1+N) % N, 2*N + i])

    for i in range(N):
        faces.append([N + i, 3*N + i, 6*N + (i + 1) % N])
        faces.append([N + i, 6*N + (i + 1) % N, 7*N + (i + 1) % N])
        faces.append([2*N + i, 7*N + i, 5*N + i])
        faces.append([2*N + i, 5*N + i, 3*N + i])

    return points, faces


def print_increments(j_no, pgon, ang, ang_end, steps, flags, **kwargs):
    p_ang = pgon.angle()/2
    ang_inc = (ang_end - ang)/steps
    for i in range(steps):
        a = ang+i*ang_inc
        try:
            if j_no == '88':
                principal_pts = j88_get_principal_verts(pgon, a, flags)
                edge_len = (principal_pts[5] - principal_pts[4]).mag()
            elif j_no == '89':
                principal_pts = j89_get_principal_verts(pgon, a, flags)
                edge_len = (principal_pts[5] - principal_pts[3]).mag()
            elif j_no == '90':
                principal_pts = j90_get_principal_verts(pgon, a, flags)
                B_rot = principal_pts[1].rot_z(-p_ang/2)
                edge_len = Vec(0, 2*B_rot[1], 2*B_rot[2]).mag()
        except:
            edge_len = -0.001
        print("%-20.16g,%.16g" % (a*180/math.pi, edge_len), **kwargs)


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'johnson_no',
        help='johnson solid number to make the general form of',
        choices=['88', '89', '90'],
        nargs='?',
        default='90')
    parser.add_argument(
        'polygon_fraction',
        help='number of sides of the base polygon (N), '
             'or a fraction for star polygons (N/D)',
        type=anti_lib.read_polygon,
        nargs='?',
        default='3')
    parser.add_argument(
        'angle',
        help='angle to rotate a belt edge from horizontal',
        type=float,
        nargs='?',
        default=0.0)
    parser.add_argument(
        'angle_end',
        help='if given, do not output the model, but instead '
             'print the free edge length for a range of models '
             'between angle and angle_end. This may be used to '
             'search for, and zero in on, unit edged models.',
        type=float,
        nargs='?',
        default=None)
    parser.add_argument(
        'steps',
        help='number of steps in the range of models when '
             'printing the free edge lengths',
        type=anti_lib.read_positive_int,
        nargs='?',
        default=20)
    parser.add_argument(
        '-p', '--print-edge',
        help='print the free edge length to the screen when the '
             'output type is model (otherwise ignore)',
        action='store_true')
    parser.add_argument(
        '-f', '--form-flags',
        help='letters from A, B, C, ... which select an '
             'alternative form at the stages of the calculation '
             'where a choice exists. These will generally "pop" '
             'a vertex type in or outfile. For each model the '
             'available letters and an example vertex popped '
             'are: 88: A-3N, B-4N; 89: A-4N, B-6N, C-5N; '
             '90: A-3N, B-N.',
        default='')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    ang = math.radians(args.angle)
    pgon = args.polygon_fraction
    flags = args.form_flags

    try:
        if args.angle_end is None:
            if args.johnson_no == '88':
                (points, faces) = j88_make_poly(pgon, ang, flags)
                edge_len = (points[5*pgon.N] - points[4*pgon.N]).mag()
            elif args.johnson_no == '89':
                (points, faces) = j89_make_poly(pgon, ang, flags)
                edge_len = (points[6*pgon.N] - points[4*pgon.N]).mag()
            elif args.johnson_no == '90':
                (points, faces) = j90_make_poly(pgon, ang, flags)
                B_rot = points[pgon.N].rot_z(-pgon.angle()/2)
                edge_len = Vec(0, 2*B_rot[1], 2*B_rot[2]).mag()
            if(args.print_edge):
                print("free edge length:", edge_len, file=sys.stderr)
            out = anti_lib.OffFile(args.outfile)
            out.print_all_pgon(points, faces, pgon)
        else:
            end_ang = math.radians(args.angle_end)
            print_increments(args.johnson_no, pgon, ang, end_ang, args.steps,
                             flags, file=args.outfile)

    except ValueError as e:
        parser.error(e.args[0])

if __name__ == "__main__":
    main()
