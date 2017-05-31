README
======

antiprism\_python
-----------------

| This is a collection of geometry modelling programs written in
  Python3,
| and associated with the `Antiprism <http://www.antiprism.com>`__
  project.
| The Antiprism programs may be used to view, process or analyse these
  models.

| The individual program scripts depend on the included anti\_lib.py
  module.
| If you are downloading a single script then make sure to also download
| **`anti\_lib.py <./anti_lib.py>`__**, and put both in the same
  directory.

| Some of these programs were written to solve a specific problem,
| some to solve a general problem, some were written as prototypes.
| The programs vary in quality, completeness and usefulness. They are
| all shared under the MIT licence.

Programs
~~~~~~~~

To see the help for a program, run it with *-h*.

-  **`barrel.py <./barrel.py>`__**
   Create a cyclic polyhedron with one or two bands of equatorial
   squares, oriented like diamonds.
-  **`eq\_antherm.py <./eq_antherm.py>`__**
   Create an antihermaphrodite with equilateral triangles
-  **`geodesic.py <./geodesic.py>`__**
   Create coordinates for a higher frequency, plane-faced or
   spherical, icosahedron, octahedron or tetrahedron.
-  **`gold\_bowtie.py <./gold_bowtie.py>`__**
   Create a polyhedron with axial symmetry involving a golden
   trapezium bow-tie
-  **`lat\_grid.py <./lat_grid.py>`__**
   Make a variety of lattices and grids using integer coordinates.
-  **`njitterbug.py <./njitterbug.py>`__**
   Create a jitterbug model for a general polygon base. The
   transformation includes, if constructible, models corresponding to
   the antiprism, snub-antiprisms and gyrobicupola for that base.
-  **`njohnson.py <./njohnson.py>`__**
   Create a Johnson-based model, from J88, J89, J90, with a general
   polygon base.
-  **`packer.py <./packer.py>`__**
   Pack balls in a sphere. The pack is seeded with two or more balls,
   then subsequent balls are added one at a time in three point contact
   in positions chosen by the packing method.
-  **`pentabelt.py <./pentabelt.py>`__**
   Make a model with axial symmetry based on a belt of staggered
   pentagons
-  **`proj\_dome.py <./proj_dome.py>`__**
   Make a Jacoby-style dome, as described in
   http://www.google.com/patents/US7900405 . Project a tiling of
   unit-edged triangles, squares or crossed squares (unit edges), at
   a specified height, onto a unit hemisphere, by gnomonic,
   stereographic
   or general point projection.
-  **`ring\_place.py <./ring_place.py>`__**
   Place maximum radius rings of contacting balls around points on a
   sphere. Input is a list of coordinates, one set per line.
-  **`sph\_circles.py <./sph_circles.py>`__**
   Distribute points on horizontal circles on a sphere (like a disco
   ball). The sphere is split into equal width bands. Balls with a
   diameter of this width are distributed equally around each band.
   The number of balls is either as many points as will fit in the band,
   or a specified number.
-  **`sph\_saff.py <./sph_saff.py>`__**
   Distribute num\_points (default 20) on a sphere using the algorithm
   from "Distributing many points on a sphere" by E.B. Saff and
   A.B.J. Kuijlaars, Mathematical Intelligencer 19.1 (1997) 5--11.
-  **`sph\_spiral.py <./sph_spiral.py>`__**
   Distribute points in a spiral on a sphere.
-  **`spiro.py <./spiro.py>`__**
   Create a spirograph pattern.
-  **`str\_art.py <./str_art.py>`__**
   Create simple epicycloid string art patterns.
-  **`temcor\_dome.py <./temcor_dome.py>`__**
   Make a Temcor-style dome, using the method described in
   https://groups.google.com/d/msg/geodesichelp/hJ3V9Nfp3kE/nikgoBPSFfwJ
   .
   The base model is a pyramid with a unit edge base polygon at a
   specified height above the origin. The axis to rotate the plane
   about passes through the origin and is in the direction of the
   base polygon mid-edge to the pyramid apex.
-  **`tri\_tiling.py <./tri_tiling.py>`__**
   Create a polyhedron which tiles the sphere with congruent triangles.
-  **`twister.py <./twister.py>`__**
   Twist two polygons placed on symmetry axes and joined by a vertex
-  **`twister\_rhomb.py <./twister_rhomb.py>`__**
   Twist polygons, of the same type, placed on certain fixed axes and
   joined by vertices.
-  **`twister\_test.py <./twister_test.py>`__**
   Twist two polygons placed on axes at a specified angle and joined by
   a vertex.

Complementary Programs
~~~~~~~~~~~~~~~~~~~~~~

Related Python programs in external projects

-  **`antitile <https://github.com/brsr/antitile>`__**
   Generates geodesic models by various methods.
-  **`view\_off.py <https://github.com/brsr/antitile/blob/master/bin/view_off.py>`__**
   An OFF file viewer with export to PNG and SVG.
