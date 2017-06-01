# README

## antiprism_python

This is a collection of geometry modelling programs written in Python,
and associated with the [Antiprism](http://www.antiprism.com) project.
The Antiprism programs may be used to view, process or analyse these models.

Some of these programs were written to solve a specific problem,
some to solve a general problem, some were written as prototypes.
The programs vary in quality, completeness and usefulness. They are
shared under the MIT licence.

### Programs

To see the help for a program, run it with *-h*.

* **[barrel.py](./barrel.py)**  
    Create a cyclic polyhedron with one or two bands of equatorial
    squares, oriented like diamonds.
* **[eq_antherm.py](./eq_antherm.py)**  
    Create an antihermaphrodite with equilateral triangles
* **[geodesic.py](./geodesic.py)**  
    Create coordinates for a higher frequency, plane-faced or
    spherical, icosahedron, octahedron or tetrahedron.
* **[gold_bowtie.py](./gold_bowtie.py)**  
    Create a polyhedron with axial symmetry involving a golden
    trapezium bow-tie
* **[lat_grid.py](./lat_grid.py)**  
    Make a variety of lattices and grids using integer coordinates.
* **[njitterbug.py](./njitterbug.py)**  
    Create a jitterbug model for a general polygon base. The
    transformation includes, if constructible, models corresponding to
    the antiprism, snub-antiprisms and gyrobicupola for that base.  
* **[njohnson.py](./njohnson.py)**  
    Create a Johnson-based model, from J88, J89, J90, with a general
    polygon base.
* **[packer.py](./packer.py)**  
    Pack balls in a sphere. The pack is seeded with two or more balls,
    then subsequent balls are added one at a time in three point contact
    in positions chosen by the packing method.
* **[pentabelt.py](./pentabelt.py)**  
    Make a model with axial symmetry based on a belt of staggered pentagons
* **[proj_dome.py](./proj_dome.py)**  
    Make a Jacoby-style dome, as described in
    http://www.google.com/patents/US7900405 . Project a tiling of
    unit-edged triangles, squares or crossed squares (unit edges), at
    a specified height, onto a unit hemisphere, by gnomonic, stereographic
    or general point projection.
* **[ring_place.py](./ring_place.py)**  
    Place maximum radius rings of contacting balls around points on a
    sphere. Input is a list of coordinates, one set per line.
* **[sph_circles.py](./sph_circles.py)**  
    Distribute points on horizontal circles on a sphere (like a disco
    ball). The sphere is split into equal width bands. Balls with a
    diameter of this width are distributed equally around each band.
    The number of balls is either as many points as will fit in the band,
    or a specified number.
* **[sph_saff.py](./sph_saff.py)**  
    Distribute num_points (default 20) on a sphere using the algorithm
    from "Distributing many points on a sphere" by E.B. Saff and
    A.B.J. Kuijlaars, Mathematical Intelligencer 19.1 (1997) 5--11.
* **[sph_spiral.py](./sph_spiral.py)**  
    Distribute points in a spiral on a sphere.
* **[spiro.py](./spiro.py)**  
    Create a spirograph pattern.
* **[str_art.py](./str_art.py)**  
    Create simple epicycloid string art patterns.
* **[temcor_dome.py](./temcor_dome.py)**  
    Make a Temcor-style dome, using the method described in
    https://groups.google.com/d/msg/geodesichelp/hJ3V9Nfp3kE/nikgoBPSFfwJ .
    The base model is a pyramid with a unit edge base polygon at a
    specified height above the origin. The axis to rotate the plane
    about passes through the origin and is in the direction of the
    base polygon mid-edge to the pyramid apex.
* **[tri_tiling.py](./tri_tiling.py)**  
    Create a polyhedron which tiles the sphere with congruent triangles.
* **[twister.py](./twister.py)**  
    Twist two polygons placed on symmetry axes and joined by a vertex
* **[twister_rhomb.py](./twister_rhomb.py)**  
    Twist polygons, of the same type, placed on certain fixed axes and
    joined by vertices.
* **[twister_test.py](./twister_test.py)**  
    Twist two polygons placed on axes at a specified angle and joined by
    a vertex.

### Complementary Programs

Related Python programs in external projects

* **[antitile](https://github.com/brsr/antitile)**  
    Generates geodesic models by various methods.
* **[view_off.py](https://github.com/brsr/antitile/blob/master/bin/view_off.py)**  
    An OFF file viewer with export to PNG and SVG.
