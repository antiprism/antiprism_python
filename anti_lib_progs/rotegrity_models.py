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
Create cyclic rotegrity models with 1 or 2 layers of units
'''
import argparse
import sys
from math import cos, sin, tan, acos, asin, pi
import anti_lib


def frac_to_end_ang(frac, num_layers):
    return (pi/num_layers) * frac / (1 - frac)


def pgon_ang_to_end_ang(a):
    return acos(2*cos(a)-1)
# return 2*asin(sqrt(1-cos(a)))
    # return 2*atan(tan(a)/sqrt(2))


def make_model(pgon, frac, num_layers):
    """Make cyclic rotegrity model"""
    N = pgon.N
    a = pgon.angle() / 2
    if num_layers == 1:
        e_ang = frac_to_end_ang(frac, num_layers)
    else:
        e_ang = pgon_ang_to_end_ang(a)

    # print("e_ang="+str(e_ang), file=sys.stderr)
    rot = asin(tan(e_ang/2)/tan(a))
    ax_ang = rot - pi/2
    mid_edge = anti_lib.Vec(sin(rot), 0, cos(rot))
    axis = anti_lib.Vec(sin(ax_ang), 0, cos(ax_ang))
    pt_angs = [-e_ang/2, e_ang/2, pi/num_layers-e_ang/2, pi/num_layers+e_ang/2]
    pts = [anti_lib.Mat.rot_axis_ang(axis, pt_ang) * mid_edge
           for pt_ang in pt_angs]
    if num_layers == 2:
        mid_ht = (pts[2][2]+pts[3][2])/2
        for i in range(4):
            pts[i][2] -= mid_ht

    points = [anti_lib.Vec(0, 0, 0)] * (2*num_layers) * N
    faces = []
    for i in range(N):
        ang = i * 2 * a
        rot_i = anti_lib.Mat.rot_xyz(0, 0, ang)
        rot_i2 = anti_lib.Mat.rot_xyz(0, 0, ang + a)
        if num_layers == 1:
            points[i] = rot_i * pts[0]
            points[i+N] = rot_i * pts[2]
            faces.append([N + (i+1) % N, i, (i+1) % N, N + i])
        else:  # num_layers == 2:
            points[i] = rot_i * pts[0]
            points[i+N] = rot_i * pts[2]
            points[i+2*N] = rot_i * pts[3]
            points[i+3*N] = rot_i2 * pts[0]
            points[i+3*N][2] *= -1

            faces.append([i+2*N, i, (i+1) % N, i+N])
            faces.append([(i+1) % N + N, i + 3*N, (i+1) % N + 3*N, i + 2*N])
            # 30,31,20,11
    return points, faces, e_ang, pi/num_layers + e_ang


def main():
    """Entry point"""
    epilog = '''
notes:
  Depends on anti_lib.py.

examples:
  octagonal hosohedron rotegrity
  rotegrity_pretwist.py 8 | antiview
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
        '-f', '--fraction',
        help='end fraction (default: 0.1), only used for 1-layer models',
        type=float,
        default='0.1')
    parser.add_argument(
        '-n', '--number-layers',
        help='number of layers (default: 1)',
        choices=['1', '2'],
        default='1')
    parser.add_argument(
        '-c', '--curved-struts',
        help='represent struts by a curved sequence of edges',
        action='store_true')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()
    number_layers = int(args.number_layers)

    pgon = args.polygon_fraction
    points, faces, end_ang, total_ang = make_model(
            pgon, args.fraction, number_layers)

    use_curves = args.curved_struts

    invisible = 9999
    out = anti_lib.OffFile(args.outfile)
    if not use_curves:
        out.print_header(len(points), 5 * len(faces))
        out.print_verts(points)
        for i in range(len(faces)):
            out.print_face(faces[i], 0, i)
            for j in range(4):
                col = i if j else invisible
                out.print_face([faces[i][j], faces[i][(j+1) % 4]], 0, col)
    else:
        num_seg_pts = 60
        out.print_header(num_seg_pts*len(faces), (num_seg_pts-1)*len(faces))
        for i in range(len(faces)):
            pts = [points[idx] for idx in faces[i]]
            axis = anti_lib.Vec.cross(pts[1]-pts[0], pts[2]-pts[0]).unit()
            for j in range(num_seg_pts):
                rot = anti_lib.Mat.rot_axis_ang(
                        axis, total_ang * j / (num_seg_pts-1))
                out.print_vert(rot * pts[1])
        for i in range(len(faces)):
            start_idx = i * num_seg_pts
            for j in range(num_seg_pts-1):
                out.print_face([start_idx+j, start_idx+j+1], 0, i)

    print("end_angle=%19.17f\ntot_angle=%19.17f\nfraction=%19.17f\n" % (
        end_ang, total_ang, end_ang/total_ang), file=sys.stderr)


if __name__ == "__main__":
    main()
