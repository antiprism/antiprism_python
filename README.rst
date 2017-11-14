README
======

antiprism_python
-----------------

`antiprism_python <https://github.com/antiprism/antiprism_python/>`_
is a collection of geometry modelling programs written in Python3,
and associated with the `Antiprism <http://www.antiprism.com>`_ project.
The Antiprism programs may be used to view, process or analyse these models.

Some of these programs were written to solve a specific problem,
some to solve a general problem, some were written as prototypes.
The programs vary in quality, completeness and usefulness. They are
all shared under the MIT licence.

Install
~~~~~~~

Install the whole package with ``pip3`` (or ``pip`` if that installs as
Python3), either from the PyPI repository::

   pip3 install antiprism_python

Or directly from the Git repository::

   pip3 install git+git://github.com/antiprism/antiprism_python

Alternatively, download the programs you are interested in and copy 
`anti_lib.py <https://github.com/antiprism/antiprism_python/blob/master/anti_lib.py>`_
into the same directory.

Programs
~~~~~~~~

To see the help for a program, run it with *-h*.

- `barrel.py <https://github.com/antiprism/antiprism_python/blob/master/barrel.py>`_
   Create a cyclic polyhedron with one or two bands of equatorial
   squares, oriented like diamonds.
- `eq_antherm.py <https://github.com/antiprism/antiprism_python/blob/master/eq_antherm.py>`_
   Create an antihermaphrodite with equilateral triangles
- `geodesic.py <https://github.com/antiprism/antiprism_python/blob/master/geodesic.py>`_
   Create coordinates for a higher frequency, plane-faced or
   spherical, icosahedron, octahedron or tetrahedron.
- `gold_bowtie.py <https://github.com/antiprism/antiprism_python/blob/master/gold_bowtie.py>`_
   Create a polyhedron with axial symmetry involving a golden
   trapezium bow-tie
- `lat_grid.py <https://github.com/antiprism/antiprism_python/blob/master/lat_grid.py>`_
   Make a variety of lattices and grids using integer coordinates.
- `lamella.py <https://github.com/antiprism/antiprism_python/blob/master/lamella.py>`_
   Create a lamella domes. Also, makes square barrel models and
   multiply-gyroelongated antihermaphrodite models.
- `njitterbug.py <https://github.com/antiprism/antiprism_python/blob/master/njitterbug.py>`_
   Create a jitterbug model for a general polygon base. The
   transformation includes, if constructible, models corresponding to
   the antiprism, snub-antiprisms and gyrobicupola for that base.
- `njohnson.py <https://github.com/antiprism/antiprism_python/blob/master/njohnson.py>`_
   Create a Johnson-based model, from J88, J89, J90, with a general
   polygon base.
- `packer.py <https://github.com/antiprism/antiprism_python/blob/master/packer.py>`_
   Pack balls in a sphere. The pack is seeded with two or more balls,
   then subsequent balls are added one at a time in three point contact
   in positions chosen by the packing method.
- `pentabelt.py <https://github.com/antiprism/antiprism_python/blob/master/pentabelt.py>`_
   Make a model with axial symmetry based on a belt of staggered
   pentagons
- `proj_dome.py <https://github.com/antiprism/antiprism_python/blob/master/proj_dome.py>`_
   Make a Jacoby-style dome, as described in
   http://www.google.com/patents/US7900405 . Project a tiling of
   unit-edged triangles, squares or crossed squares (unit edges), at
   a specified height, onto a unit hemisphere, by gnomonic,
   stereographic
   or general point projection.
- `ring_place.py <https://github.com/antiprism/antiprism_python/blob/master/ring_place.py>`_
   Place maximum radius rings of contacting balls around points on a
   sphere. Input is a list of coordinates, one set per line.
- `sph_circles.py <https://github.com/antiprism/antiprism_python/blob/master/sph_circles.py>`_
   Distribute points on horizontal circles on a sphere (like a disco
   ball). The sphere is split into equal width bands. Balls with a
   diameter of this width are distributed equally around each band.
   The number of balls is either as many points as will fit in the band,
   or a specified number.
- `sph_saff.py <https://github.com/antiprism/antiprism_python/blob/master/sph_saff.py>`_
   Distribute num_points (default 20) on a sphere using the algorithm
   from "Distributing many points on a sphere" by E.B. Saff and
   A.B.J. Kuijlaars, Mathematical Intelligencer 19.1 (1997) 5--11.
- `sph_spiral.py <https://github.com/antiprism/antiprism_python/blob/master/sph_spiral.py>`_
   Distribute points in a spiral on a sphere.
- `spiro.py <https://github.com/antiprism/antiprism_python/blob/master/spiro.py>`_
   Create a spirograph pattern.
- `str_art.py <https://github.com/antiprism/antiprism_python/blob/master/str_art.py>`_
   Create simple epicycloid string art patterns.
- `temcor_dome.py <https://github.com/antiprism/antiprism_python/blob/master/temcor_dome.py>`_
   Make a Temcor-style dome, using the method described in
   https://groups.google.com/d/msg/geodesichelp/hJ3V9Nfp3kE/nikgoBPSFfwJ
   .
   The base model is a pyramid with a unit edge base polygon at a
   specified height above the origin. The axis to rotate the plane
   about passes through the origin and is in the direction of the
   base polygon mid-edge to the pyramid apex.
- `tri_tiling.py <https://github.com/antiprism/antiprism_python/blob/master/tri_tiling.py>`_
   Create a polyhedron which tiles the sphere with congruent triangles.
- `twister.py <https://github.com/antiprism/antiprism_python/blob/master/twister.py>`_
   Twist two polygons placed on symmetry axes and joined by a vertex
- `twister_rhomb.py <https://github.com/antiprism/antiprism_python/blob/master/twister_rhomb.py>`_
   Twist polygons, of the same type, placed on certain fixed axes and
   joined by vertices.
- `twister_test.py <https://github.com/antiprism/antiprism_python/blob/master/twister_test.py>`_
   Twist two polygons placed on axes at a specified angle and joined by
   a vertex.

Complementary Programs
~~~~~~~~~~~~~~~~~~~~~~

Related Python programs in external projects

- `antitile <https://github.com/brsr/antitile>`_
   Generates geodesic models by various methods.
- `view_off.py <https://github.com/brsr/antitile/blob/master/bin/view_off.py>`_
   An OFF file viewer with export to PNG and SVG.
