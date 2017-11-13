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
anti_lib.py: private library of common code used by programs in this
repository.
'''


import argparse
import fractions
import math
import sys


def safe_for_trig(val):
    if abs(val) > 1:
        return -1 if val < 0 else 1
    else:
        return val


class Vec:
    def __init__(self, *v):
        self.v = list(v)

    def fromlist(self, v):
        if not isinstance(v, list):
            raise TypeError
        self.v = v[:]
        return self

    def copy(self):
        return Vec().fromlist(self.v)

    def __str__(self):
        return '(' + repr(self.v)[1:-1] + ')'

    def __repr__(self):
        return 'Vec(' + repr(self.v)[1:-1] + ')'

    def __len__(self):
        return len(self.v)

    def __getitem__(self, key):
        if not isinstance(key, int):
            raise TypeError
        if key < 0 or key >= len(self.v):
            raise KeyError
        return self.v[key]

    def __setitem__(self, key, value):
        if not isinstance(key, int):
            raise TypeError
        if key < 0 or key >= len(self.v):
            raise KeyError
        self.v[key] = value

    # Element-wise negation
    def __neg__(self):
        v = list(map(lambda x: -x, self.v))
        return Vec().fromlist(v)

    # Element-wise addition
    def __add__(self, other):
        v = list(map(lambda x, y: x+y, self.v, other.v))
        return Vec().fromlist(v)

    # Element-wise subtraction
    def __sub__(self, other):
        v = list(map(lambda x, y: x-y, self.v, other.v))
        return Vec().fromlist(v)

    # Element-wise multiplication by scalar
    def __mul__(self, scalar):
        v = list(map(lambda x: x*scalar, self.v))
        return Vec().fromlist(v)

    # Element-wise pre-multiplication by scalar
    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    # Element-wise division by scalar
    def __truediv__(self, scalar):
        return self.__mul__(1/scalar)

    # Vector magnitude/length squared
    def mag2(self):
        return self.dot(self, self)

    # Vector magnitude/length
    def mag(self):
        return math.sqrt(self.mag2())

    # Vector as unit
    def unit(self):
        return self.__truediv__(self.mag())

    # Vector rotated about z-axis
    def rot_z(self, ang):
        r = math.sqrt(self.v[0]**2+self.v[1]**2)
        initial_ang = math.atan2(self.v[1], self.v[0])
        final_ang = initial_ang + ang
        return Vec(r*math.cos(final_ang), r*math.sin(final_ang), self.v[2])

    # Cross product v0 x v1
    @staticmethod
    def cross(v0, v1):
        return Vec(v1[2]*v0[1] - v1[1]*v0[2],
                   v1[0]*v0[2] - v1[2]*v0[0],
                   v1[1]*v0[0] - v1[0]*v0[1])

    # Dot product v0 . v1
    @staticmethod
    def dot(v0, v1):
        return sum(map(lambda x, y: x*y, v0.v, v1.v))

    # Triple product v0. (v1 x v2)
    @staticmethod
    def triple(v0, v1, v2):
        return Vec.dot(v0, Vec.cross(v1, v2))


def centroid(v):
    return sum(v, Vec(0, 0, 0)) / len(v)


class Mat:
    def __init__(self, init_type=1):
        if init_type == 0:
            self.m = [0]*16
        elif init_type == 1:
            self.m = [0]*16
            for i in range(4):
                self.m[5*i] = 1

    @staticmethod
    def rot_axis_ang(axis, ang):
        rot = Mat()
        c = math.cos(ang)
        s = math.sin(ang)
        axis = axis.unit()
        t = 1.0 - c
        rot.m[0] = c + axis[0]*axis[0]*t
        rot.m[5] = c + axis[1]*axis[1]*t
        rot.m[10] = c + axis[2]*axis[2]*t

        tmp1 = axis[0]*axis[1]*t
        tmp2 = axis[2]*s
        rot.m[4] = tmp1 + tmp2
        rot.m[1] = tmp1 - tmp2

        tmp1 = axis[0]*axis[2]*t
        tmp2 = axis[1]*s
        rot.m[8] = tmp1 - tmp2
        rot.m[2] = tmp1 + tmp2

        tmp1 = axis[1]*axis[2]*t
        tmp2 = axis[0]*s
        rot.m[9] = tmp1 + tmp2
        rot.m[6] = tmp1 - tmp2

        return rot

    @staticmethod
    def rot_xyz(ang_x, ang_y, ang_z):
        mat = Mat()
        if ang_z:
            mat *= Mat.rot_axis_ang(Vec(0, 0, 1), ang_z)
        if ang_y:
            mat *= Mat.rot_axis_ang(Vec(0, 1, 0), ang_y)
        if ang_x:
            mat *= Mat.rot_axis_ang(Vec(1, 0, 0), ang_x)

        return mat

    @staticmethod
    def rot_from_to(v_from, v_to):
        epsilon = 1e-13
        n_from = v_from.unit()
        n_to = v_to.unit()
        axis = Vec.cross(n_from, n_to)
        cos_a = safe_for_trig(Vec.dot(n_from, n_to))
        if abs(cos_a) >= 1-epsilon:
            axis = Vec.cross(n_from, Vec(1.2135, 2.09865, 3.23784)).unit()

        return Mat.rot_axis_ang(axis, math.acos(cos_a))

    @staticmethod
    def rot_from_to2(from1, from2, to1, to2):
        r = (Mat.rot_xyz(0.1, 0, 0) * Mat.rot_xyz(0, 0.2, 0) *
             Mat.rot_xyz(0, 0, 0.3))
        inv_r = (Mat.rot_xyz(0, 0, -0.3) * Mat.rot_xyz(0, -0.2, 0) *
                 Mat.rot_xyz(-0.1, 0, 0))
        to1 = r * to1
        to2 = r * to2

        trans = Mat.rot_from_to(from1, to1)
        from1 = trans * from1
        from2 = trans * from2

        norm1 = Vec.cross(from2, from1).unit()
        norm2 = Vec.cross(to2, to1).unit()

        # find the angle to rotate about to1, in the range -180<ang<=180
        ang = math.acos(safe_for_trig(Vec.dot(norm1, norm2)))
        if Vec.triple(to1, norm1, norm2) < 0:
            ang *= -1

        trans = Mat.rot_axis_ang(to1, ang) * trans
        return inv_r*trans

    @staticmethod
    def transl(vec):
        mat = Mat()
        mat.m[3], mat.m[7], mat.m[11] = vec[0], vec[1], vec[2]
        return mat

    # Multiplication matrix * vector
    def mult_vec(self, vec):
        new_point = Vec(0, 0, 0)
        for i in range(12):
            if i % 4 != 3:
                new_point[i//4] += self.m[i] * vec[i % 4]
            else:
                new_point[i//4] += self.m[i]

        return new_point

    # Multiplication matrix * matrix
    def mult_matrix(self, other):
        ret = Mat(0)
        for i in range(16):
            for j in range(4):
                ret.m[i] += self.m[(i//4)*4+j]*other.m[j*4 + (i % 4)]
        return ret

    # Multiplication matrix * (matrix or vector)
    def __mul__(self, other):
        if isinstance(other, Mat):
            return self.mult_matrix(other)
        elif isinstance(other, Vec):
            return self.mult_vec(other)
        else:
            raise Exception('Matrix multiplication with unsupported type')


def angle_around_axis(v0, v1, axis):
    n0 = Vec.cross(v0, axis).unit()
    n1 = Vec.cross(v1, axis).unit()
    ang = math.acos(safe_for_trig(Vec.dot(n0, n1)))
    if Vec.dot(axis, Vec.cross(n0, n1)) < 0:
        ang = 2*math.pi - ang
    return ang


class RawFraction:
    def __init__(self, N, D=None):
        if D is None:
            D = 1
        elif D == 0:
            raise ValueError('fraction denominator is zero')
        if N < 0:
            N *= -1
            D *= -1
        self.N = N
        self.D = D

    def read(self, frac):
        D = 1
        vals = frac.split('/')
        sz = len(vals)
        if sz > 2:
            raise ValueError('fraction has more than one /')
        if sz > 1:
            try:
                D = int(vals[1])
            except:
                raise ValueError('fraction denominator is not an integer')
        try:
            N = int(vals[0])
        except:
            raise ValueError('fraction numerator is not an integer')

        self.__init__(N, D)


class Polygon(RawFraction):
    def __init__(self, N=2, D=None):
        if D is None:
            D = 1
        if N < 0:
            N *= -1
            D *= -1
        if N < 2:
            raise ValueError('fraction numerator is less than 2')

        D %= N
        if D == 0:
            raise ValueError('fraction denominator is a multiple '
                             'of the numerator')
        self.parts = fractions.gcd(N, D)
        self.N = N // self.parts
        self.D = D // self.parts

    def angle(self):
        return 2*math.pi*self.D/self.N

    def circumradius(self, edge=1.0):
        return edge / (2*math.sin(self.angle()/2))

    def inradius(self, edge=1.0):
        return edge / (2*math.tan(self.angle()/2))

    def get_points(self, P=Vec(1, 0, 0)):
        return [P.rot_z((j*self.angle() + i*2*math.pi/(self.N*self.parts)))
                for i in range(self.parts) for j in range(self.N)]

    def get_faces(self, offset=0):
        return [[i+j*self.N+offset for i in range(self.N)]
                for j in range(self.parts)]


class OffFile:
    def __init__(self, strm=sys.stdout):
        self.strm = strm

    def print_header(self, num_verts, num_faces):
        print('OFF\n{} {} 0'.format(num_verts, num_faces), file=self.strm)

    def print_vert(self, vert):
        print(vert[0], vert[1], vert[2], file=self.strm)

    def print_face(self, face, offset=0, col=None):
        print(len(face), *[v+offset for v in face], end='', file=self.strm)
        if type(col) is int:
            print('', col, file=self.strm)
        else:
            print(file=self.strm)

    def print_verts(self, verts):
        for v in verts:
            self.print_vert(v)

    def print_faces(self, faces, offset=0, col=None):
        for face in faces:
            self.print_face(face, offset, col)

    def print_all(self, verts, faces):
        self.print_header(len(verts), len(faces))
        self.print_verts(verts)
        self.print_faces(faces)

    def print_all_pgon(self, verts, faces, pgon, repeat_side=False):
        sides = repeat_side and pgon.N or 1
        parts = pgon.parts
        self.print_header(len(verts)*sides*parts, len(faces)*sides*parts)
        for i in range(parts):
            for side in range(sides):
                trans = Mat.rot_axis_ang(Vec(0, 0, 1), 2*math.pi*(
                    side/sides + (i/parts)/pgon.N))
                self.print_verts([trans * v for v in verts])

        for i in range(parts*sides):
            self.print_faces(faces, len(verts)*i)


def read_positive_float(str_val):
    try:
        val = float(str_val)
    except:
        raise argparse.ArgumentTypeError('not a number')

    if val <= 0.0:
        raise argparse.ArgumentTypeError('not a positive number')

    return val


def read_positive_int(str_val, min_val=1):
    try:
        val = int(str_val)
    except:
        raise argparse.ArgumentTypeError('not an integer')

    if val < min_val:
        msgs = ['not zero or greater', 'not a positive integer']
        raise argparse.ArgumentTypeError(msgs[min_val])

    return val


def read_polygon(val_str):
    pgon = Polygon()
    try:
        pgon.read(val_str)
    except ValueError as e:
        raise argparse.ArgumentTypeError(e.args[0])
    return pgon


class DefFormatter(argparse.RawDescriptionHelpFormatter):
    '''Allow individual options to be pre-formatted. Description and
       epilog are unformatted'''
    def _split_lines(self, text, width):
        if text.startswith(']'):
            return text[1:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)
