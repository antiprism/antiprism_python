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
Make a lattice or grid using integer coordinates. NOTE: if this program
does not meet all your needs then also consider the Antiprism lat_grid
and lat_util programs http://www.antiprism.com/programs/
'''

import sys
import argparse
import anti_lib
from anti_lib import Vec

epsilon = 1e-8


def dist2_inside(v1, width):
    """Square of distance less than"""
    return (v1.mag2() < width+epsilon)


def dist2_test(v1, v2, idx1, idx2, len2):
    """Square of distance equal"""
    return (v1-v2).mag2() == len2


def sc_coord_test(x, y, z):  # dist2 = 1, 2, 3
    """Test for coordinate in SC lattice"""
    return 1


def fcc_coord_test(x, y, z):  # dist2 = 2, 4, 6, 12
    """Test for coordinate in FCC lattice"""
    return (x+y+z) % 2 == 0


def bcc_coord_test(x, y, z):  # dist2 = 3, 4, 8
    """Test for coordinate in BCC lattice"""
    return (x % 2 == y % 2) and (y % 2 == z % 2)


def rh_dodec_coord_test(x, y, z):  # dist2 = 3 (8)
    """Test for coordinate in rhombic cuboctahedron grid"""
    return ((x % 2 + y % 2 + z % 2) == 0 or
            ((x % 2 + y % 2 + z % 2) == 3 and (x//2 + y//2) % 2 == (z//2) % 2))


def cubo_oct_coord_test(x, y, z):  # dist2 = 2
    """Test for coordinate in octahedron/cuboctahedron grid"""
    return (x % 2 + y % 2 + z % 2) == 2


def tr_oct_coord_test(x, y, z):  # dist2 = 2
    """Test for coordinate in truncated octahedron grid"""
    return ((z % 4 == 0 and x % 4 and y % 4 and (x-y) % 2) or
            (y % 4 == 0 and z % 4 and x % 4 and (z-x) % 2) or
            (x % 4 == 0 and y % 4 and z % 4 and (y-z) % 2))


def tr_tet_tet_coord_test(x, y, z):  # dist2 = 2
    """Test for coordinate in tretrahedron/truncated tetrahedron grid"""
    return ((x % 2 == 0 and y % 4 == ((x + z)) % 4) or
            (y % 2 == 0 and z % 4 == ((y + x)) % 4) or
            (z % 2 == 0 and x % 4 == ((z + y)) % 4))


def tr_tet_tr_oct_cubo_coord_test(x, y, z):  # dist2 = 4
    """Test for coordinate in truncated tetrahedron/truncated octahedron/
    cuboctahedron grid"""
    x = abs(x) % 6
    y = abs(y) % 6
    z = abs(z) % 6
    if x > 3:
        x = 6-x
    if y > 3:
        y = 6-y
    if z > 3:
        z = 6-z
    dist2 = x**2 + y**2
    return ((z % 6 == 0 and (dist2 == 2 or dist2 == 8)) or
            (z % 6 == 1 and (dist2 == 1 or dist2 == 13)) or
            (z % 6 == 2 and (dist2 == 4 or dist2 == 10)) or
            (z % 6 == 3 and dist2 == 5))


def diamond_coord_test(x, y, z):  # dist2 = 3
    """Test for coordinate in diamond grid"""
    return (((x % 2 + y % 2 + z % 2) == 0 and (x//2+y//2+z//2) % 2 == 0) or
            ((x % 2 + y % 2 + z % 2) == 3 and (x//2+y//2+z//2) % 2 == 0))


def hcp2_coord_test(x, y, z):
    """Test for coordinate in HCP alternative grid"""
    return (x+y+z) % 6 == 0 and ((x-y) % 3 or (x-z) % 3)


def hcp_coord_test(x, y, z):
    """Test for coordinate in HCP grid"""
    m = x+y+z
    return m % 6 == 0 and (x-m//12) % 3 == 0 and (y-m//12) % 3 == 0


# Coordinates from Wendy Krieger
def hcp_diamond_coord_test(x, y, z):
    """Test for coordinate in HCP diamond grid"""
    for pt in [[0, 0, 0], [3, 3, 3], [6, 0, 6], [9, 3, 9]]:
        tx = x - pt[0]
        ty = y - pt[1]
        tz = z - pt[2]
        s = tx + ty + tz
        if(s % 24 == 0):
            n8 = s//3
            if (tx-n8) % 6 == 0 and (ty-n8) % 6 == 0 and (tz-n8) % 6 == 0:
                return True
    return False


def make_verts(lat_type, width, container_is_cube):
    """Generate coordinates"""
    verts = []
    if container_is_cube:
        coords_start = 0
    else:  # sphere
        coords_start = -width
    coords_end = width

    coord_test = eval(lat_type+"_coord_test")
    for x in range(coords_start, coords_end+1):
        for y in range(coords_start, coords_end+1):
            for z in range(coords_start, coords_end+1):
                if container_is_cube or dist2_inside(Vec(x, y, z), width):
                    if coord_test(x, y, z):
                        verts.append(Vec(x, y, z))
    return verts


def make_edges(verts, strut_len_sqd):
    """Add edges"""
    edges = []
    if strut_len_sqd:
        for i in range(len(verts)-1):
            for j in range(i+1, len(verts)):
                if dist2_test(verts[i], verts[j], i, j, strut_len_sqd):
                    edges.append((i, j))
    return edges


def main():
    """Entry point"""
    epilog = '''
notes:
  Depends on anti_lib.py.

examples:
  Oct-tet truss, space-frame
  lat_grid.py fcc 5 2 | antiview

  Spherical section of diamond
  lat_grid.py -s diamond 20 3 | antiview

'''

    parser = argparse.ArgumentParser(formatter_class=anti_lib.DefFormatter,
                                     description=__doc__, epilog=epilog)

    parser.add_argument(
        'lat_type',
        help=''']Type of lattice or grid (default: sc). The numbers in
brackets are suitable strut arguments.

sc                 - Simple Cubic             (1, 2, 3)
fcc                - Face Centred Cubic   (2, 4, 6, 12)
bcc                - Body Centred Cubic       (3, 4, 8)
hcp                - Hexagonal Close Packing       (18)
rh_dodec           - Rhombic Dodecahedra         (3, 8)
cubo_oct           - Cuboctahedron / Octahedron     (2)
tr_oct             - Truncated Octahedron           (2)
tr_tet_tet         - Trunc Tetrahedron/Tetrahedron  (2)
tr_tet_tr_oct_cubo - Trunc Tet / Trunc Oct / Cuboct (4)
diamond            - Diamond                        (3)
hcp_diamond        - HCP Diamond                   (27)''',
        choices=['sc', 'fcc', 'bcc', 'hcp', 'rh_dodec', 'cubo_oct',
                 'tr_oct', 'tr_tet_tet', 'tr_tet_tr_oct_cubo',
                 'diamond', 'hcp_diamond'],
        nargs='?',
        default="sc")
    parser.add_argument(
        'width',
        help='''width of container cube (default: 6), or radius
of contianer sphere''',
        type=anti_lib.read_positive_int,
        nargs='?',
        default=6)
    parser.add_argument(
        'strut_len_sqd',
        help='''add struts between pairs of vertices that are
separated by the square root of this distance''',
        type=anti_lib.read_positive_int,
        nargs='?',
        default=0)
    parser.add_argument(
        '-s', '--sphere-container',
        help='''container is a sphere (radius 'width') rather than
a cube''',
        action='store_true')
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    verts = make_verts(args.lat_type, args.width, not args.sphere_container)
    edges = make_edges(verts, args.strut_len_sqd)

    out = anti_lib.OffFile(args.outfile)
    out.print_all(verts, edges)

if __name__ == "__main__":
    main()
