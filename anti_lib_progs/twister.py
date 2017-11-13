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
Twist two polygons placed on symmetry axes and joined by a vertex. Output
model in OFF format.
'''

import argparse
import sys
import math
import re
import anti_lib
from anti_lib import Vec, Mat

epsilon = 1e-13


def make_arc(v0, v1, num_segs, start_idx):
    axis = Vec.cross(v0, v1).unit()
    ang = anti_lib.angle_around_axis(v0, v1, axis)
    points = [v0]
    faces = []
    mat = Mat.rot_axis_ang(axis, ang/num_segs)
    for i in range(num_segs):
        # accumulated error
        points.append(mat * points[-1])
        faces.append([start_idx + i, start_idx + i+1])
    return points, faces


def make_frame(frame_elems, axis_angle, rad, num_segs):
    points = []
    faces = []
    if frame_elems:
        v0 = Vec(0, 0, -1)
        v1 = Vec(math.sin(axis_angle), 0, -math.cos(axis_angle))
        if 'r' in frame_elems:
            ps, fs = make_arc(v0, v1, num_segs, 0)
            points += ps
            faces += fs
        if 'a' in frame_elems:
            num_pts = len(points)
            points += [v0, -v0, v1, -v1]
            faces += [[num_pts+0, num_pts+1], [num_pts+2, num_pts+3]]

        points = [rad * p for p in points]
        faces += [[i] for i in range(len(points))]

    return points, faces


def parse_axis_pair(axes_str):
    pat = re.compile('(?P<sym_type>[^[]+)'
                     '\[(?P<axis1>[^,]+),(?P<axis2>[^\]]+)](?P<id_no>.*)')
    match = pat.match(axes_str)
    if not match:
        raise argparse.ArgumentTypeError(
            'axis pair is not in form sym_type[axis1,axis2]id_no '
            'e.g I[5,3], D5[2,2]2')

    axis_pair = {}
    sym_type = match.group('sym_type')
    axis_pair['sym_type_str'] = sym_type   # for reporting
    if sym_type in ['I', 'O', 'T']:
        axis_pair['sym_type'] = sym_type
    elif sym_type[0] == "D":
        axis_pair['sym_type'] = 'D'
        if len(sym_type) == 1:
            raise argparse.ArgumentTypeError(
                'n-fold number missing for D symmetry')
        try:
            axis_pair['dih_N'] = int(sym_type[1:])
        except:
            raise argparse.ArgumentTypeError(
                'n-fold number for D symmetry is not an integer')
        if axis_pair['dih_N'] < 2:
            raise argparse.ArgumentTypeError(
                'n-fold number for D symmetry cannot be less than 2')
    else:
        raise argparse.ArgumentTypeError(
            'invalid symmetry type \'%s\', should be from T, O, I, or can be '
            'D followed by an integer or fraction (e.g D5) '
            % (sym_type))

    axis_pair['nfolds'] = []
    for i, axis in enumerate(['axis1', 'axis2']):
        try:
            axis_pair['nfolds'].append(int(match.group(axis)))
        except ValueError:
            raise argparse.ArgumentTypeError(
                ['first', 'second'][i] + ' axis is \'' + match.group(axis) +
                '\', should be a positive integer ')
        if axis_pair['nfolds'][i] < 2:
            raise argparse.ArgumentTypeError(
                ['first', 'second'][i] + ' axis cannot be less than 2')

    if not match.group('id_no'):
        axis_pair['id_no'] = 1
    else:
        try:
            axis_pair['id_no'] = int(match.group('id_no'))
        except ValueError:
            raise argparse.ArgumentTypeError(
                'id_no is \'' + match.group('id_no') +
                '\', should be a positive integer ')
        if axis_pair['id_no'] < 1:
            raise argparse.ArgumentTypeError('id_no cannot be less than 1')

    return axis_pair


def read_axes(axes_str):
    axes_str = ''.join(axes_str.split())   # remove whitespace
    axis_pair = parse_axis_pair(axes_str)
    # place greater nfold axis first
    if axis_pair['nfolds'][0] >= axis_pair['nfolds'][1]:
        nfolds = tuple(axis_pair['nfolds'])
        axes_in_order = True
    else:
        nfolds = tuple(axis_pair['nfolds'][::-1])
        axes_in_order = False

    phi = (math.sqrt(5)+1)/2
    A = {}
    if axis_pair['sym_type'] == 'T':
        # T: (3, 3), (3, 2), (2, 2)
        A[(3, 3)] = [
            [[1, 1, 1], [-1, -1, 1]]]
        A[(3, 2)] = [
            [[1, 1, 1], [0, 0, 1]]]
        A[(2, 2)] = [
            [[0, 0, 1], [1, 0, 0]]]

    elif axis_pair['sym_type'] == 'O':
        # O: (4, 4), (4, 3), (4, 2)*2, (3, 3), (3, 2)*2, (2, 2)*2
        A[(4, 4)] = [
            [[0, 0, 1], [1, 0, 0]]]                 # 1
        A[(4, 3)] = [
            [[0, 0, 1], [1, 1, 1]]]                 # 1
        A[(4, 2)] = [
            [[0, 0, 1], [0, 1, 1]],                 # 1
            [[0, 0, 1], [1, 1, 0]]]                 # 2
        A[(3, 3)] = [
            [[1, 1, 1], [1, -1, 1]]]                # 1
        A[(3, 2)] = [
            [[1, 1, 1], [0, -1, -1]],               # 1
            [[1, 1, 1], [1, 0, -1]]]                # 2
        A[(2, 2)] = [
            [[0, 1, 1], [1, 0, 1]],                 # 1
            [[0, 1, 1], [0, 1, -1]]]                # 2

    elif axis_pair['sym_type'] == 'I':
        # I: (5, 5), (5, 3)*2, (5, 2)*3, (3, 3)*2, (3, 2)*4, (2, 2)*4
        A[(5, 5)] = [
            [[0, 1, phi], [0, 1, -phi]]]            # 1
        A[(5, 3)] = [
            [[0, 1, phi], [1, 1, 1]],               # 1
            [[0, 1, phi], [phi, -1/phi, 0]]]        # 2
        A[(5, 2)] = [
            [[0, 1, phi], [0, 0, -1]],              # 1
            [[0, 1, phi], [1, 1/phi, -phi]],        # 2
            [[0, 1, phi], [1, 0, 0]]]               # 3
        A[(3, 3)] = [
            [[1, 1, 1], [-1/phi, 0, -phi]],         # 1
            [[1, 1, 1], [1, -1, -1]]]               # 2
        A[(3, 2)] = [
            [[1, 1, 1], [-1, -1/phi, -phi]],        # 1
            [[1, 1, 1], [-1, 0, 0]],                # 2
            [[1, 1, 1], [1, -1/phi, -phi]],         # 3
            [[1, 1, 1], [1, 1/phi, -phi]]]          # 4
        A[(2, 2)] = [
            [[0, 0, 1], [1, 1/phi, phi]],           # 1
            [[0, 0, 1], [1/phi, phi, 1]],           # 2
            [[0, 0, 1], [phi, 1, 1/phi]],           # 3
            [[0, 0, 1], [1, 0, 0]]]                 # 4

    else:  # axis_pair['sym_type'] == 'D':
        N = axis_pair['dih_N']
        A[(N, 2)] = [[[0, 0, 1], [1, 0, 0]]]
        A[(2, 2)] = [[[1, 0, 0],
                     [math.cos(i*math.pi/N), math.sin(i*math.pi/N), 0]]
                     for i in range(1, 1 + N//2)]

    if nfolds not in A:
        raise argparse.ArgumentTypeError(
            '[%d, %d] is not a valid axis pair for symmetry type %s' % (
                axis_pair['nfolds'][0], axis_pair['nfolds'][1],
                axis_pair['sym_type_str']))

    if axis_pair['id_no'] > len(A[nfolds]):
        raise argparse.ArgumentTypeError(
            'id_no for symmetry axes %s[%d, %d] must be less than %d' % (
                axis_pair['sym_type_str'],
                axis_pair['nfolds'][0], axis_pair['nfolds'][1],
                len(A[nfolds])))

    axis_pair['axes'] = [
        Vec().fromlist(axis).unit() for axis in
        A[nfolds][axis_pair['id_no']-1]]
    if not axes_in_order:
        axis_pair['axes'] = [-axis_pair['axes'][1], -axis_pair['axes'][0]]

    return axis_pair


def calc_polygons(pgon0, pgon1, ang, ratio, angle_between_axes):
    # rotate initial vertex of polygon0 by ang. Common
    # point lies on a line through vertex in the direction
    # of axis0 (z-axis). Common vertex lies on a cylinder
    # radius r1 with axis on axis1.

    # rotate model so axis1 is on z-axis, to simplify the
    # calculation of the intersection
    rot = Mat.rot_axis_ang(Vec(0, 1, 0), angle_between_axes)

    # initial vertex turned on axis
    V = Vec(pgon0.circumradius(), 0, 0).rot_z(ang)
    # vertex in rotated model
    Q = rot * V
    # direction of line in rotated model
    u = rot * Vec(0, 0, 1)

    # equation of line is P = Q + tu for components x y z and parameter t
    # equation of cylinder at z=0 is x^2 + y^2 = r1^2
    a = u[0]**2 + u[1]**2
    b = 2*(Q[0]*u[0] + Q[1]*u[1])
    c = Q[0]**2 + Q[1]**2 - (pgon1.circumradius()*ratio)**2

    disc = b**2 - 4*a*c
    if disc < -epsilon:
        raise Exception("model is not geometrically constructible")
    elif disc < 0:
        disc = 0

    t = (-b - math.sqrt(disc))/(2*a)

    P = V + Vec(0, 0, t)   # the common point

    points = pgon0.get_points(P)
    faces = pgon0.get_faces()

    Q = rot * P
    rot_inv = Mat.rot_axis_ang(Vec(0, 1, 0), -angle_between_axes)
    points += [rot_inv*p for p in pgon1.get_points(Q)]
    faces += pgon1.get_faces(pgon0.N*pgon0.parts)

    return points, faces


def read_axis_multiplier(mult_str):
    mult_fract = anti_lib.RawFraction(1)
    if len(mult_str) > 0:
        try:
            mult_fract.read(mult_str)
        except ValueError as e:
            raise argparse.ArgumentTypeError(
                'invalid multiplier fraction: ' + e.args[0])

    return mult_fract


def frame_type(arg):
    if arg.strip('ra'):
        raise argparse.ArgumentTypeError(
            'frame type contains letters other than r, a')
    return arg


def read_turn_angle(ang_str):
    if ang_str[-1] == 'e':
        ang_type = 'e'
        ang_str = ang_str[:-1]
    elif ang_str[-1] == 'x':
        ang_type = 'x'
        ang_str = ang_str[:-1]
    else:
        ang_type = 'd'
    try:
        ang_val = float(ang_str)
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            'invalid numeric value: ', e.args[0])

    return [ang_val, ang_type]


def main():
    """Entry point"""

    epilog = '''
notes:
  Depends on anti_lib.py. Use poly_kscope to repeat the model.

examples:
  Icosahedral model
  twister.py I[5,3] | poly_kscope -s I| antiview

  Twisted icosahedral model with hexagons
  twister.py I[5,3] 1 2 -a 0.5e | poly_kscope -s I| antiview

  Dihedral model with frame
  twister.py D6[6,2] 1 2 -f ra | poly_kscope -s D6 | antiview
'''

    parser = argparse.ArgumentParser(formatter_class=anti_lib.DefFormatter,
                                     description=__doc__, epilog=epilog)

    parser.add_argument(
        'symmetry_axes',
        help=']Axes given in the form: sym_type[axis1,axis2]id_no\n'
             '   sym_type: rotational symmetry group, can be T, O, I,\n'
             '       or can be D followed by an integer (e.g D5)\n'
             '   axis1,axis2: rotational order of each of the two axes\n'
             '   id_no (default: 1): integer to select between\n'
             '       non-equivalent pairs of axes having the same\n'
             '       symmetry group and rotational orders\n'
             'e.g. T[2,3], I[5,2]2, D7[7,2], D11[2,2]4\n'
             'Axis pairs are from the following\n'
             '   T: [3, 3], [3, 2], [2, 2]\n'
             '   O: [4, 4], [4, 3], [4, 2]x2, [3, 3], [3, 2]x2,\n'
             '      [2, 2]x2\n'
             '   I: [5, 5], [5, 3]x2, [5, 2]x3, [3, 3]x2, [3, 2]x4,\n'
             '      [2, 2]x4\n'
             '   Dn: [n 2], [2,2]x(n/2 rounded down)\n',
        type=read_axes,
        nargs='?',
        default='O[4,3]')
    parser.add_argument(
        'multiplier1',
        help='integer or fractional multiplier for axis 1 '
             '(default: 1). If the axis is N/D and the '
             'multiplier is n/d the polygon used will be N*n/d',
        type=read_axis_multiplier,
        nargs='?',
        default="1")
    parser.add_argument(
        'multiplier2',
        help='integer or fractional multiplier for axis 2 '
             '(default: 1). If the axis is N/D and the '
             'multiplier is n/d the polygon used will be N*n/d',
        type=read_axis_multiplier,
        nargs='?',
        default="1")
    parser.add_argument(
        '-r', '--ratio',
        help='ratio of edge lengths (default: 1.0)',
        type=float,
        default=1.0)
    parser.add_argument(
        '-a', '--angle',
        help='amount to turn polygon on axis0 in degrees '
        '(default: 0), or a value followed by \'e\', '
             'where 1.0e is half the central angle of an edge, '
             'which produces an edge-connected model (negative '
             'values may have to be specified as, e.g. -a=-1.0e), '
             'or a value followed by x, which is like e but with '
             'a half turn offset',
        type=read_turn_angle,
        default='0')
    parser.add_argument(
        '-f', '--offset',
        help='amount to offset the first polygon to avoid '
             'coplanarity with the second polygon, for example '
             '0.0001 (default: 0.0)',
        type=float,
        default=0.0)
    parser.add_argument(
        '-F', '--frame',
        help='include frame elements in output, any from: '
             'r - rhombic tiling edges, '
             'a - rotation axes (default: no elements)',
        type=frame_type,
        default='')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    axis_pair = args.symmetry_axes
    pgons = []
    for i, m in enumerate([args.multiplier1, args.multiplier2]):
        try:
            pgons.append(anti_lib.Polygon(
                axis_pair['nfolds'][i]*m.N, m.D))
        except Exception as e:
                parser.error('multiplier%d: ' % (i) + e.args[0])

    axes = axis_pair['axes']

    if args.angle[1] == 'e':     # units: half edge central angle
        turn_angle = args.angle[0] * pgons[0].angle()/2
    elif args.angle[1] == 'x':   # units: half edge central angle
        turn_angle = math.pi + args.angle[0] * pgons[0].angle()/2
    else:                        # units: degrees
        turn_angle = math.radians(args.angle[0])

    sin_angle_between_axes = Vec.cross(axes[0], axes[1]).mag()
    if abs(sin_angle_between_axes) > 1:
        sin_angle_between_axes = -1 if sin_angle_between_axes < 0 else 1
    angle_between_axes = math.asin(sin_angle_between_axes)
    if(Vec.dot(axes[0], axes[1]) > 0):
        axes[1] *= -1

    try:
        (points, faces) = calc_polygons(
            pgons[0], pgons[1], turn_angle, args.ratio, angle_between_axes)
    except Exception as e:
        parser.error(e.args[0])

    if args.offset:
        for i in range(len(faces[0])):
            points[i][2] += args.offset

    frame_rad = calc_polygons(pgons[0], pgons[1], 0, args.ratio,
                              angle_between_axes)[0][0].mag()
    frame_points, frame_faces = make_frame(args.frame, angle_between_axes,
                                           frame_rad, 10)

    rot = Mat.rot_from_to2(Vec(0, 0, 1), Vec(1, 0, 0), axes[0], axes[1])
    all_points = [rot * point for point in points+frame_points]

    out = anti_lib.OffFile(args.outfile)
    out.print_header(len(all_points), len(faces)+len(frame_faces))
    out.print_verts(all_points)
    for i in range(pgons[0].parts+pgons[1].parts):
        out.print_face(faces[i], 0, int(i < pgons[0].parts))
    out.print_faces(frame_faces, len(points), 2)


if __name__ == "__main__":
    main()
