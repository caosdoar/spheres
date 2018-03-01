# 4 Ways to Create a Mesh for a Sphere

## Introduction

When graphics programmers face the problem of creating a mesh for a sphere, trade-offs must be made between quality and construction, memory and rendering costs. This document introduces four different methods and analyzes their characteristics and compares them allowing programmers to make an informed decision on which method suits their needs.

## Building the mesh

![alt text](https://github.com/caosdoar/spheres/raw/master/img/spheres_4.png "Spheres comparison")

### Standard Sphere

This is the most common implementation of a sphere mesh and can be found in almost any 3d toolset. You can find it under the name of “UV sphere” in blender or just “sphere” in 3d max. This method divides the sphere using meridians (lines from pole to pole) and parallels (lines parallel to the equator). It produces faces with bigger area near the equator and smaller ones close to the poles.

```
for j in parallels_count:
  parallel = PI * (j+1) / parallels_count
  for i in meridians_count:
    meridian = 2.0 * PI * i / meridians_count
    return spherical_to_cartesian(meridian, parallel)
```

### Normalized Cube

This method uses a uniformly subdivided cube where each vertex position is normalized and multiplied by the sphere radius. This creates a non uniformly subdivided sphere where the triangles closer to the center of a cube face are bigger than the ones closer to the edges of the cube.

```
for f in faces:
	  origin = get_origin(f)
	  right = get_right_dir(f)
	  up = get_up_dir(f)
	  for j in subdiv_count:
	    for i in subdiv_count:
	      face_point = origin + 2.0 * (right * i + up * j) / subdiv_count
	      return normalize(face_point)
```

### Spherified Cube

This method is based on a subdivided cube as well but it tries to create more uniform divisions in the sphere. The area of the faces and the length of the edges suffer less variation, but the sphere still has some obvious deformation as points get closer to the corners of the original cube.

```
for f in faces:
  origin = get_origin(f)
  right = get_right_dir(f)
  up = get_up_dir(f)
  for j in subdiv_count:
    for i in subdiv_count:
      p = origin + 2.0 * (right * i + up * j) / subdiv_count
      p2 = p * p
      rx = sqrt(1.0 - 0.5 * (p2.y + p2.z) + p2.y*p2.z/3.0)
      ry = sqrt(1.0 - 0.5 * (p2.z + p2.x) + p2.z*p2.x/3.0)
      rz = sqrt(1.0 - 0.5 * (p2.x + p2.y) + p2.x*p2.y/3.0)
      return (rx, ry, rz)
```

### Icosahedron

An icosahedron is a polyhedron composed of 20 identical equilateral triangles and possesses some interesting properties: Each triangle has the same area and each vertex is at the same distance from all its neighbours. 

To get a higher number of triangles we need to subdivide each triangle into four triangles by creating a new vertex at the middle point of each edge which is then normalized, to make it lie in the sphere surface. Sadly this breaks the initial properties of the icosahedron, the triangles are not equilateral anymore and neither the area nor the distance between adjacent vertices is the same across the mesh. An added problem with this method is that we can only increase the number of faces by four each time. But it is still a better approximation by almost any measure excluding its number of triangles.

The pseudocode is not added here due to the length of the algorithm that needs the initial 12 vertices and 20 faces to be manually written. The pseudocode for the subdivision algorithm looks like the following:

```
for f in input_mesh.faces:
  v0 = input_mesh.vertex(f, 0)
  v1 = input_mesh.vertex(f, 1)
  v2 = input_mesh.vertex(f, 2)
  v3 = normalize(0.5 * (v0 + v1))
  v4 = normalize(0.5 * (v1 + v2))
  v5 = normalize(0.5 * (v2 + v0))
  output_mesh.add_face(v0, v3, v5)
  output_mesh.add_face(v3, v1, v4)
  output_mesh.add_face(v5, v4, v2)
  output_mesh.add_face(v3, v4, v5)
```

## How is a Mesh Better than Another

The mesh will be an imperfect representation as we are using triangles to approximate the surface of the sphere. The first metric we will use to compare the different implementations is the distance from a number of points on the sphere surface to the mesh created.

Another thing we are going to test is the ratio between the expected triangle area and the actual triangle area. This can be interesting for subdivision purposes as more uniform triangles may give better results.

Both metrics will be given as maximum and mean square error. Smaller numbers are better.

## Results

All charts in this section show the number of triangles in the X axis and the (average or maximum) error in the Y axis.

### Distance to the Sphere: Average Error

![alt text](https://github.com/caosdoar/spheres/raw/master/img/surface_avg_error.png "Surface average error")

As we can see the standard sphere has the largest average error and the icosahedron is better than any of the other methods but for the number of triangles analyzed here we can only create two subdivisions with 320 and 1280 triangles, which limits the flexibility of the method.

### Distance to the Sphere: Maximum Error

![alt text](https://github.com/caosdoar/spheres/raw/master/img/surface_max_error.png "Surface maximum error")

Despite the average error to be quite close between the normalized and the spherified cube we can see here that the maximum error is larger for the first one been worst than in the uv sphere. The icosahedron still wins in accuracy over the other three.

### Area per Triangle: Mean Square Error

![alt text](https://github.com/caosdoar/spheres/raw/master/img/area_avg_error.png "Surface area average error")

Here is where we see the benefits of the spherified cube. It is the only one where the ratio between the triangles area and the expected triangle area (area of the sphere divided by the number of triangles in the mesh) decreases with the number of triangles.The error in the icosahedron case increases because successive subdivisions produces for each triangle four triangles of slightly different areas.

### Area per Triangle: Maximum Square Error

![alt text](https://github.com/caosdoar/spheres/raw/master/img/area_max_error.png "Surface area maximum error")

We can see that the behaviour is similar to the maximum error of surface distances. When looking at the average error, the normalized cube performs better than the UV sphere but the situation is reversed when considering the maximum error. In this case the icosahedron beats the spherified cube for all the analyzed subdivisions but we can see how the error is already increasing.

## Visualization of the surface distance error

![alt text](https://github.com/caosdoar/spheres/raw/master/img/errors_l.jpg "Surface errors in the spheres")

## Notes 

The implementation provided for the spheres based on subdivided cubes generates duplicated vertices on the edges of the cube. This could be avoided at the cost of complicating the code.

The triangulation of the subdivided cubes can be improved by making them symmetrical to the center of the faces of the cube.

![alt text](https://github.com/caosdoar/spheres/raw/master/img/correct_triangulation.png "Correct way of triangulating a spherical cube")

## Conclusion

The UV sphere is the worst option in terms of accuracy for a given number of triangles but it is still the easiest algorithm. If you are going to render a lot of spheres or you need highly tessellated spheres think about switching to any of the other algorithms.

## Thanks

Thanks to @ppalmier for supervising this document.
