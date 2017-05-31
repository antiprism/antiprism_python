#!/usr/bin/env python3

# Copyright (c) 2003-2016 Adrian Rossiter <adrian@antiprism.com>
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
Distribute points on horizontal circles on a sphere (like a disco ball).
The sphere is split into equal width bands. Balls with a diameter of this
width are distributed equally around each band. The number of balls is
either as many points as will fit in the band, or a specified number.
The output is a list of ball centres.
'''

import argparse
import sys
import math
import anti_lib
from anti_lib import Vec


def calc_points(args):
    points = []
    vert_ang_inc = math.pi/(args.number_circles-1)
    ball_dia = math.sqrt(2 - 2*math.cos(vert_ang_inc))
    num_circles = args.number_circles
    horz_stagger = 0.0
    circles = []
    for circle in range(num_circles):
        if args.points_list:
            num_balls = args.points_list[circle]
            if num_balls == -1:
                num_balls = args.default_number_points
        else:
            num_balls = args.default_number_points

        vert_ang = circle*vert_ang_inc
        if num_balls == -1:
            rad = math.sin(vert_ang)
            try:
                num_balls = int(math.floor(2*math.pi /
                                math.acos((2*rad*rad-ball_dia*ball_dia) /
                                          (2*rad*rad)) * args.aspect_ratio))
            except:
                num_balls = 1

        horz_ang_inc = 2*math.pi/num_balls if num_balls > 0 else 0
        circles.append([vert_ang, num_balls, horz_ang_inc, 0.0])

    if args.stagger:
        for circle in range(num_circles//2, -1, -1):
            if circle == num_circles//2:
                circles[circle][3] = circles[circle][2]/4
                circles[circle+1][3] = -circles[circle][2]/4
            else:
                for circ in [[circle, circle+1],
                             [num_circles-circle-1, num_circles-circle-2]]:
                    horz_stagger = circles[circ[1]][3]
                    if horz_stagger > 0:
                        horz_stagger -= circles[circ[1]][2]/2
                    else:
                        horz_stagger += circles[circ[1]][2]/2
                    circles[circ[0]][3] = horz_stagger

    for circle in range(args.exclude_poles, num_circles-args.exclude_poles):
        vert_ang = circles[circle][0]
        num_balls = circles[circle][1]
        horz_ang_inc = circles[circle][2]
        horz_stagger = circles[circle][3]
        rad = math.sin(vert_ang)

        for n in range(num_balls):
            horz_ang = n*horz_ang_inc+horz_stagger
            points.append(Vec(rad*math.cos(horz_ang), rad*math.sin(horz_ang),
                              math.cos(vert_ang)))

    return points


def read_nonnegative_int_list(str_val):
    val = []
    for cnt in str_val.split(','):
        if cnt == '':
            val.append(-1)
        else:
            try:
                val.append(anti_lib.read_positive_int(cnt, min_val=0))
            except:
                raise argparse.ArgumentTypeError(
                    'circle point counts must be non-negative integers '
                    'separated only by commas')
    return val


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(description=__doc__)

    num_group = parser.add_mutually_exclusive_group()
    num_group.add_argument(
        '-n', '--number-circles',
        help='number of circles of points (default: 6)',
        type=anti_lib.read_positive_int,
        default=6)
    num_group.add_argument(
        '-p', '--points-list',
        help='number of points on each circles',
        type=read_nonnegative_int_list)
    parser.add_argument(
        '-d', '--default-number-points',
        help='number of points to put on a circle if not '
             'otherwise specified (default: maximum number '
             'that fit)',
        type=anti_lib.read_positive_int,
        default=-1)
    parser.add_argument(
        '-s', '--stagger',
        help='stagger points between layers',
        action='store_true')
    parser.add_argument(
        '-x', '--exclude-poles',
        help='exclude the pole point circles',
        action='store_true')
    parser.add_argument(
        '-a', '--aspect-ratio',
        help='angular aspect ratio (default: 1.0)',
        type=anti_lib.read_positive_float,
        default=1.0)
    parser.add_argument(
        '-o', '--outfile',
        help='output file name (default: standard output)',
        type=argparse.FileType('w'),
        default=sys.stdout)

    args = parser.parse_args()

    if args.points_list:
        args.number_circles = len(args.points_list)

    points = calc_points(args)

    out = anti_lib.OffFile(args.outfile)
    out.print_all(points, [])

if __name__ == "__main__":
    main()
