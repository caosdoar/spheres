// Oscar Sebio Cajaraville 2015
// caosdoar@gmail.com
//
// What sphere is the best?
// Or...
// For a number of triangles what sphere is closest to the perfect sphere?

/*
MIT License

Copyright (c) 2015 Oscar Sebio Cajaraville

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#define _CRT_SECURE_NO_WARNINGS

#define EXPORT_OBJS 0
#define EXPORT_IMAGES 0
#define PRINT_STATS 1
#define PRINT_CREATION_TIMES 0

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <map>
#include <random>
#include <string>
#include <iostream>
#include <fstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#if _WIN32
#include <Windows.h>
#else
typedef uint64_t LARGE_INTEGER;
#endif

struct Vector3
{
	double x, y, z;
	Vector3(double v) : x(v), y(v), z(v) {}
	Vector3(double ix, double iy, double iz) : x(ix), y(iy), z(iz) {}
	Vector3 operator +(const Vector3 &other) const { return Vector3(x + other.x, y + other.y, z + other.z); }
	Vector3 operator -(const Vector3 &other) const { return Vector3(x - other.x, y - other.y, z - other.z); }
	Vector3 operator *(const Vector3 &other) const { return Vector3(x * other.x, y * other.y, z * other.z); }
};

double dot(const Vector3 &a, const Vector3 &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 cross(const Vector3 &a, const Vector3 &b)
{
	return Vector3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

double length(const Vector3 &a)
{
	return std::sqrt(dot(a, a));
}

Vector3 normalize(const Vector3 &a)
{
	const double lrcp = 1.0 / std::sqrt(dot(a, a));
	return Vector3(a.x * lrcp, a.y * lrcp, a.z * lrcp);
}

struct Mesh
{
	std::vector<Vector3> vertices;
	std::vector<uint32_t> triangles;

	uint32_t triangleCount() const { return triangles.size() / 3; }

	void addTriangle(uint32_t a, uint32_t b, uint32_t c)
	{
		triangles.emplace_back(a);
		triangles.emplace_back(b);
		triangles.emplace_back(c);
	}

	void addQuad(uint32_t a, uint32_t b, uint32_t c, uint32_t d)
	{
		triangles.emplace_back(a);
		triangles.emplace_back(b);
		triangles.emplace_back(c);
		triangles.emplace_back(a);
		triangles.emplace_back(c);
		triangles.emplace_back(d);
	}

	void addQuadAlt(uint32_t a, uint32_t b, uint32_t c, uint32_t d)
	{
		triangles.emplace_back(a);
		triangles.emplace_back(b);
		triangles.emplace_back(d);
		triangles.emplace_back(b);
		triangles.emplace_back(c);
		triangles.emplace_back(d);
	}

	void clear()
	{
		vertices.clear();
		triangles.clear();
	}

	double distance(const Vector3 &p, uint32_t tidx) const
	{
		const uint32_t idx0 = triangles[tidx];
		const uint32_t idx1 = triangles[tidx + 1];
		const uint32_t idx2 = triangles[tidx + 2];
		const Vector3 v0 = vertices[idx0];
		const Vector3 v1 = vertices[idx1];
		const Vector3 v2 = vertices[idx2];
		const Vector3 bv = v0;
		const Vector3 e0 = v1 - v0;
		const Vector3 e1 = v2 - v0;
		const Vector3 dv = bv - p;
		const double a = dot(e0, e0);
		const double b = dot(e0, e1);
		const double c = dot(e1, e1);
		const double d = dot(e0, dv);
		const double e = dot(e1, dv);
		const double f = dot(dv, dv);

		const double det = a*c - b*b;
		double s = b*e - c*d;
		double t = b*d - a*e;

		if (s + t <= det)
		{
			if (s < 0.0)
			{
				if (t < 0.0)
				{
					// region 4
					if (d < 0.0)
					{
						t = 0.0;
						s = -d >= a ? 1.0 : -d / a;
					}
					else
					{
						s = 0.0;
						t = e >= 0.0 ? 0.0 : (-e >= c ? 1.0 : -e / c);
					}
				}
				else
				{
					// region 3
					s = 0.0;
					t = e >= 0.0 ? 0.0 : (-e >= c ? 1.0 : -e / c);
				}
			}
			else if (t < 0.0)
			{
				// region 5
				s = d >= 0.0 ? 0.0 : (-d >= a ? 1.0 : -d / a);
				t = 0.0;
			}
			else
			{
				// region 0
				const double invDet = 1.0 / det;
				s *= invDet;
				t *= invDet;
			}
		}
		else
		{
			if (s < 0.0)
			{
				// region 2
				const double tmp0 = b + d;
				const double tmp1 = c + e;
				if (tmp1 > tmp0)
				{
					const double numer = tmp1 - tmp0;
					const double denom = a - 2.0 * b + c;
					s = numer >= denom ? 1.0 : numer / denom;
					t = 1.0 - s;
				}
				else
				{
					s = 0.0;
					t = (tmp1 <= 0.0 ? 1.0 : (e >= 0.0 ? 0.0 : -e / c));
				}
			}
			else if (t < 0.0)
			{
				// region 6
				const double tmp0 = b + e;
				const double tmp1 = a + d;
				if (tmp1 > tmp0)
				{
					const double numer = tmp1 - tmp0;
					const double denom = a - 2.0 * b + c;
					t = numer >= denom ? 1.0 : numer / denom;
					s = 1.0 - t;
				}
				else
				{
					s = (tmp1 <= 0.0 ? 1.0 : (d >= 0.0 ? 0.0 : -d / a));
					t = 0.0;
				}
			}
			else
			{
				// region 1
				const double numer = c + e - b - d;
				if (numer <= 0)
				{
					s = 0.0;
				}
				else
				{
					const double denom = a - 2.0 * b + c;
					s = numer >= denom ? 1.0 : numer / denom;
				}
				t = 1.0 - s;
			}
		}

		return length(p - (v0 + Vector3(s) * e0 + Vector3(t) * e1));
	}

	double distance(const Vector3 &p) const
	{
		double min = 10e10;
		for (uint32_t i = 0; i < triangles.size(); i += 3)
		{
			min = std::fmin(min, distance(p, i));
		}
		return min;
	}
};

void UVSphere(uint32_t meridians, uint32_t parallels, Mesh &mesh)
{
	mesh.vertices.emplace_back(0.0f, 1.0f, 0.0f);
	for (uint32_t j = 0; j < parallels - 1; ++j)
	{
		double const polar = M_PI * double(j+1) / double(parallels);
		double const sp = std::sin(polar);
		double const cp = std::cos(polar);
		for (uint32_t i = 0; i < meridians; ++i)
		{
			double const azimuth = 2.0 * M_PI * double(i) / double(meridians);
			double const sa = std::sin(azimuth);
			double const ca = std::cos(azimuth);
			double const x = sp * ca;
			double const y = cp;
			double const z = sp * sa;
			mesh.vertices.emplace_back(x, y, z);
		}
	}
	mesh.vertices.emplace_back(0.0f, -1.0f, 0.0f);

	for (uint32_t i = 0; i < meridians; ++i)
	{
		uint32_t const a = i + 1;
		uint32_t const b = (i + 1) % meridians + 1;
		mesh.addTriangle(0, b, a);
	}

	for (uint32_t j = 0; j < parallels - 2; ++j)
	{
		uint32_t aStart = j * meridians + 1;
		uint32_t bStart = (j + 1) * meridians + 1;
		for (uint32_t i = 0; i < meridians; ++i)
		{
			const uint32_t a = aStart + i;
			const uint32_t a1 = aStart + (i + 1) % meridians;
			const uint32_t b = bStart + i;
			const uint32_t b1 = bStart + (i + 1) % meridians;
			mesh.addQuad(a, a1, b1, b);
		}
	}

	for (uint32_t i = 0; i < meridians; ++i)
	{
		uint32_t const a = i + meridians * (parallels - 2) + 1;
		uint32_t const b = (i + 1) % meridians + meridians * (parallels - 2) + 1;
		mesh.addTriangle(mesh.vertices.size() - 1, a, b);
	}
}

namespace CubeToSphere
{
	static const Vector3 origins[6] =
	{
		Vector3(-1.0, -1.0, -1.0),
		Vector3(1.0, -1.0, -1.0),
		Vector3(1.0, -1.0, 1.0),
		Vector3(-1.0, -1.0, 1.0),
		Vector3(-1.0, 1.0, -1.0),
		Vector3(-1.0, -1.0, 1.0)
	};
	static const Vector3 rights[6] =
	{
		Vector3(2.0, 0.0, 0.0),
		Vector3(0.0, 0.0, 2.0),
		Vector3(-2.0, 0.0, 0.0),
		Vector3(0.0, 0.0, -2.0),
		Vector3(2.0, 0.0, 0.0),
		Vector3(2.0, 0.0, 0.0)
	};
	static const Vector3 ups[6] =
	{
		Vector3(0.0, 2.0, 0.0),
		Vector3(0.0, 2.0, 0.0),
		Vector3(0.0, 2.0, 0.0),
		Vector3(0.0, 2.0, 0.0),
		Vector3(0.0, 0.0, 2.0),
		Vector3(0.0, 0.0, -2.0)
	};
};

void NormalizedCube(uint32_t divisions, Mesh &mesh)
{
	const double step = 1.0 / double(divisions);
	const Vector3 step3(step, step, step);

	for (uint32_t face = 0; face < 6; ++face)
	{
		const Vector3 origin = CubeToSphere::origins[face];
		const Vector3 right = CubeToSphere::rights[face];
		const Vector3 up = CubeToSphere::ups[face];
		for (uint32_t j = 0; j < divisions + 1; ++j)
		{
			const Vector3 j3(j, j, j);
			for (uint32_t i = 0; i < divisions + 1; ++i)
			{
				const Vector3 i3(i, i, i);
				const Vector3 p = origin + step3 * (i3 * right + j3 * up);
				mesh.vertices.emplace_back(normalize(p));
			}
		}
	}

	const uint32_t k = divisions + 1;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t j = 0; j < divisions; ++j)
		{
			const bool bottom = j < (divisions / 2);
			for (uint32_t i = 0; i < divisions; ++i)
			{
				const bool left = i < (divisions / 2);
				const uint32_t a = (face * k + j) * k + i;
				const uint32_t b = (face * k + j) * k + i + 1;
				const uint32_t c = (face * k + j + 1) * k + i;
				const uint32_t d = (face * k + j + 1) * k + i + 1;
				if (bottom ^ left) mesh.addQuadAlt(a, c, d, b);
				else mesh.addQuad(a, c, d, b);
			}
		}
	}
}

void SpherifiedCube(uint32_t divisions, Mesh &mesh)
{
	const double step = 1.0 / double(divisions);
	const Vector3 step3(step, step, step);

	for (uint32_t face = 0; face < 6; ++face)
	{
		const Vector3 origin = CubeToSphere::origins[face];
		const Vector3 right = CubeToSphere::rights[face];
		const Vector3 up = CubeToSphere::ups[face];
		for (uint32_t j = 0; j < divisions + 1; ++j)
		{
			const Vector3 j3(j, j, j);
			for (uint32_t i = 0; i < divisions + 1; ++i)
			{
				const Vector3 i3(i, i, i);
				const Vector3 p = origin + step3 * (i3 * right + j3 * up);
				const Vector3 p2 = p * p;
				const Vector3 n
				(
					p.x * std::sqrt(1.0 - 0.5 * (p2.y + p2.z) + p2.y*p2.z / 3.0),
					p.y * std::sqrt(1.0 - 0.5 * (p2.z + p2.x) + p2.z*p2.x / 3.0),
					p.z * std::sqrt(1.0 - 0.5 * (p2.x + p2.y) + p2.x*p2.y / 3.0)
				);
				mesh.vertices.emplace_back(n);
			}
		}
	}

	const uint32_t k = divisions + 1;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t j = 0; j < divisions; ++j)
		{
			const bool bottom = j < (divisions / 2);
			for (uint32_t i = 0; i < divisions; ++i)
			{
				const bool left = i < (divisions / 2);
				const uint32_t a = (face * k + j) * k + i;
				const uint32_t b = (face * k + j) * k + i + 1;
				const uint32_t c = (face * k + j + 1) * k + i;
				const uint32_t d = (face * k + j + 1) * k + i + 1;
				if (bottom ^ left) mesh.addQuadAlt(a, c, d, b);
				else mesh.addQuad(a, c, d, b);
			}
		}
	}
}

void Icosahedron(Mesh &mesh)
{
	const double t = (1.0 + std::sqrt(5.0)) / 2.0;

	// Vertices
	mesh.vertices.emplace_back(normalize(Vector3(-1.0,  t, 0.0)));
	mesh.vertices.emplace_back(normalize(Vector3( 1.0,  t, 0.0)));
	mesh.vertices.emplace_back(normalize(Vector3(-1.0, -t, 0.0)));
	mesh.vertices.emplace_back(normalize(Vector3( 1.0, -t, 0.0)));
	mesh.vertices.emplace_back(normalize(Vector3(0.0, -1.0,  t)));
	mesh.vertices.emplace_back(normalize(Vector3(0.0,  1.0,  t)));
	mesh.vertices.emplace_back(normalize(Vector3(0.0, -1.0, -t)));
	mesh.vertices.emplace_back(normalize(Vector3(0.0,  1.0, -t)));
	mesh.vertices.emplace_back(normalize(Vector3( t, 0.0, -1.0)));
	mesh.vertices.emplace_back(normalize(Vector3( t, 0.0,  1.0)));
	mesh.vertices.emplace_back(normalize(Vector3(-t, 0.0, -1.0)));
	mesh.vertices.emplace_back(normalize(Vector3(-t, 0.0,  1.0)));

	// Faces
	mesh.addTriangle(0, 11, 5);
	mesh.addTriangle(0, 5, 1);
	mesh.addTriangle(0, 1, 7);
	mesh.addTriangle(0, 7, 10);
	mesh.addTriangle(0, 10, 11);
	mesh.addTriangle(1, 5, 9);
	mesh.addTriangle(5, 11, 4);
	mesh.addTriangle(11, 10, 2);
	mesh.addTriangle(10, 7, 6);
	mesh.addTriangle(7, 1, 8);
	mesh.addTriangle(3, 9, 4);
	mesh.addTriangle(3, 4, 2);
	mesh.addTriangle(3, 2, 6);
	mesh.addTriangle(3, 6, 8);
	mesh.addTriangle(3, 8, 9);
	mesh.addTriangle(4, 9, 5);
	mesh.addTriangle(2, 4, 11);
	mesh.addTriangle(6, 2, 10);
	mesh.addTriangle(8, 6, 7);
	mesh.addTriangle(9, 8, 1);
}

struct Edge
{
	uint32_t v0;
	uint32_t v1;

	Edge(uint32_t v0, uint32_t v1)
		: v0(v0 < v1 ? v0 : v1)
		, v1(v0 < v1 ? v1 : v0)
	{
	}

	bool operator <(const Edge &rhs) const
	{
		return v0 < rhs.v0 || (v0 == rhs.v0 && v1 < rhs.v1);
	}
};

uint32_t subdivideEdge(uint32_t f0, uint32_t f1, const Vector3 &v0, const Vector3 &v1, Mesh &io_mesh, std::map<Edge, uint32_t> &io_divisions)
{
	const Edge edge(f0, f1);
	auto it = io_divisions.find(edge);
	if (it != io_divisions.end())
	{
		return it->second;
	}

	const Vector3 v = normalize(Vector3(0.5) * (v0 + v1));
	const uint32_t f = io_mesh.vertices.size();
	io_mesh.vertices.emplace_back(v);
	io_divisions.emplace(edge, f);
	return f;
}

void SubdivideMesh(const Mesh &meshIn, Mesh &meshOut)
{
	meshOut.vertices = meshIn.vertices;

	std::map<Edge, uint32_t> divisions; // Edge -> new vertex

	for (uint32_t i = 0; i < meshIn.triangleCount(); ++i)
	{
		const uint32_t f0 = meshIn.triangles[i * 3];
		const uint32_t f1 = meshIn.triangles[i * 3 + 1];
		const uint32_t f2 = meshIn.triangles[i * 3 + 2];

		const Vector3 v0 = meshIn.vertices[f0];
		const Vector3 v1 = meshIn.vertices[f1];
		const Vector3 v2 = meshIn.vertices[f2];

		const uint32_t f3 = subdivideEdge(f0, f1, v0, v1, meshOut, divisions);
		const uint32_t f4 = subdivideEdge(f1, f2, v1, v2, meshOut, divisions);
		const uint32_t f5 = subdivideEdge(f2, f0, v2, v0, meshOut, divisions);

		meshOut.addTriangle(f0, f3, f5);
		meshOut.addTriangle(f3, f1, f4);
		meshOut.addTriangle(f4, f2, f5);
		meshOut.addTriangle(f3, f4, f5);
	}
}

struct Error
{
	double max;
	double mse;
};

Error surfaceError(const Mesh &mesh, uint32_t count)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1.0, 1.0);

	Error error = { 0.0, 0.0 };
	for (uint32_t i = 0; i < count; ++i)
	{
		Vector3 const p(dis(gen), dis(gen), dis(gen));
		Vector3 const pn = normalize(p);
		double e = mesh.distance(pn);
		error.max = std::fmax(e, error.max);
		error.mse += e;
	}
	error.mse /= double(count);

	return error;
}

Error surfaceAreaError(const Mesh &mesh)
{
	const double expectedTriangleArea = 4.0 * M_PI / double(mesh.triangleCount());

	Error error = { 0.0, 0.0 };
	for (uint32_t i = 0; i < mesh.triangleCount(); ++i)
	{
		const Vector3& a = mesh.vertices[mesh.triangles[i * 3 + 0]];
		const Vector3& b = mesh.vertices[mesh.triangles[i * 3 + 1]];
		const Vector3& c = mesh.vertices[mesh.triangles[i * 3 + 2]];
		const double triangleArea = 0.5 * length(cross(b - a, c - a));
		const double e = (triangleArea - expectedTriangleArea) / expectedTriangleArea;
		error.max = std::fmax(e*e, error.max);
		error.mse += e*e;
	}
	error.mse /= double(mesh.triangleCount());

	return error;
}

void exportObj(const Mesh &mesh, const std::string &name)
{
	std::fstream fs;
	fs.open(name.c_str(), std::fstream::out | std::fstream::trunc);

	for (const auto &v : mesh.vertices)
	{
		fs << "v " << v.x << " " << v.y << " " << v.z << std::endl;
	}

	for (uint32_t i = 0; i < mesh.triangleCount(); ++i)
	{
		fs << "f " << (mesh.triangles[i * 3] + 1) 
			<< " " << (mesh.triangles[i * 3 + 1] + 1)
			<< " " << (mesh.triangles[i * 3 + 2] + 1) << std::endl;
	}

	fs.close();
}

struct Pixel
{
	uint8_t r, g, b, a;
	Pixel() : r(0), g(0), b(0), a(255) {}
	Pixel(uint8_t ri, uint8_t gi, uint8_t bi, uint8_t ai = 255) : r(ri), g(gi), b(bi), a(ai) {}
};

uint8_t lerp(uint8_t a, uint8_t b, double t)
{
	const double t1 = 1.0 - t;
	return uint8_t(double(b) * t + double(a) * t1);
}

Pixel lerp(const Pixel &a, const Pixel &b, double t)
{
	return Pixel
	(
		lerp(a.r, b.r, t),
		lerp(a.g, b.g, t),
		lerp(a.b, b.b, t),
		lerp(a.a, b.a, t)
	);
}

void exportSurfaceErrorImage(const Mesh &mesh, uint32_t pointsCount, const std::string &filename, uint32_t w, uint32_t h, double maxError)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1.0, 1.0);

	uint32_t pixelsCount = w * h;
	Pixel *pixels = new Pixel[pixelsCount];
	memset(pixels, 0, sizeof(Pixel) * pixelsCount);

	const Pixel goodPixel(0, 255, 0);
	const Pixel badPixel(255, 0, 0);

	const uint32_t hw = w / 2;
	const uint32_t qw = hw / 2;
	const uint32_t hh = h / 2;
	const uint32_t cx = qw;
	const uint32_t cy = hh;
	const uint32_t cx1 = hw + qw;

	//uint32_t errorCount = 0;

	for (uint32_t i = 0; i < pointsCount; ++i)
	{
		const Vector3 p(dis(gen), dis(gen), dis(gen));
		const Vector3 pn = normalize(p);
		const double e = mesh.distance(pn);
		const double t = std::fmin(e / maxError, 1.0);
		const Pixel pixel = lerp(goodPixel, badPixel, t);

		if (pn.z > 0.0)
		{
			const uint32_t x = uint32_t(int32_t(double(qw) * pn.x) + int32_t(cx));
			const uint32_t y = uint32_t(int32_t(double(hh) * pn.y) + int32_t(cy));
			const uint32_t idx = x + y * w;
			pixels[idx] = pixel;
		}

		if (pn.y > 0.0)
		{
			const uint32_t x = uint32_t(int32_t(double(qw) * pn.x) + int32_t(cx1));
			const uint32_t y = uint32_t(int32_t(double(hh) * pn.z) + int32_t(cy));
			const uint32_t idx = x + y * w;
			pixels[idx] = pixel;
		}

		//if (e > maxError) ++errorCount;
	}

	stbi_write_png(filename.c_str(), w, h, 4, pixels, w * sizeof(Pixel));
	delete[] pixels;

	//std::cout << errorCount << std::endl;
}

void ExportMeshes()
{
	// Example on how to export it meshes in wavefront obj format
	Mesh m0;
	UVSphere(22, 11, m0);
	exportObj(m0, "uvsphere.obj");

	Mesh m1;
	NormalizedCube(6, m1);
	exportObj(m1, "ncube.obj");

	Mesh m2;
	SpherifiedCube(6, m2);
	exportObj(m2, "scube.obj");

	{
		Mesh m3;
		Mesh m4;
		Mesh m5;
		Icosahedron(m3);
		SubdivideMesh(m3, m4);
		SubdivideMesh(m4, m5);
		exportObj(m5, "icosahedron.obj");
	}
}

void ExportErrorImages()
{
	const uint32_t w = 1024;
	const uint32_t h = 512;
	const double maxError = 0.02;
	const double pointCount = 1000000;
	{
		Mesh m0;
		UVSphere(22, 11, m0);
		exportSurfaceErrorImage(m0, pointCount, std::string("uvsphere_error.png"), w, h, maxError);
	}
	{
		Mesh m0;
		NormalizedCube(6, m0);
		exportSurfaceErrorImage(m0, pointCount, std::string("normalizedcube_error.png"), w, h, maxError);
	}
	{
		Mesh m0;
		SpherifiedCube(6, m0);
		exportSurfaceErrorImage(m0, pointCount, std::string("spherifiedcube_error.png"), w, h, maxError);
	}
	{
		Mesh m3;
		Mesh m4;
		Mesh m5;
		Icosahedron(m3);
		SubdivideMesh(m3, m4);
		SubdivideMesh(m4, m5);
		exportSurfaceErrorImage(m5, pointCount, std::string("icosahedron_subdiv2_error.png"), w, h, maxError);
	}
}

void PrintStats()
{
	std::cout << "count, surface MSE, surface max. error, area per triangle MSE, area per triangle max. error" << std::endl;

	std::cout << "UVSphere" << std::endl;
	for (uint32_t i = 0; i < 22; ++i)
	{
		Mesh m;
		UVSphere(2 * (i + 2), i + 2, m);
		Error se = surfaceError(m, 10000);
		Error sae = surfaceAreaError(m);
		std::cout << m.triangleCount() << ", " << se.mse << ", " << se.max << ", " << sae.mse << ", " << sae.max << std::endl;
	}

	std::cout << "NormalizedCube" << std::endl;
	for (uint32_t i = 1; i < 14; ++i)
	{
		Mesh m;
		NormalizedCube(i, m);
		Error se = surfaceError(m, 10000);
		Error sae = surfaceAreaError(m);
		std::cout << m.triangleCount() << ", " << se.mse << ", " << se.max << ", " << sae.mse << ", " << sae.max << std::endl;
	}

	std::cout << "SpherifiedCube" << std::endl;
	for (uint32_t i = 1; i < 14; ++i)
	{
		Mesh m;
		SpherifiedCube(i, m);
		Error se = surfaceError(m, 10000);
		Error sae = surfaceAreaError(m);
		std::cout << m.triangleCount() << ", " << se.mse << ", " << se.max << ", " << sae.mse << ", " << sae.max << std::endl;
	}

	std::cout << "Icosahedron" << std::endl;
	{
		Mesh meshes[2];
		uint32_t idx = 0;
		Icosahedron(meshes[0]);
		for (uint32_t i = 0; i < 4; ++i)
		{
			Mesh &m = meshes[idx];

			Error se = surfaceError(m, 10000);
			Error sae = surfaceAreaError(m);

			std::cout << m.triangleCount() << ", " << se.mse << ", " << se.max << ", " << sae.mse << ", " << sae.max << std::endl;

			idx = (idx + 1) % 2;
			Mesh &mOut = meshes[idx];
			mOut.vertices.clear();
			mOut.triangles.clear();
			SubdivideMesh(m, mOut);
		}
	}
}

class Timer
{
public:
	void start()
	{
#if _WIN32
		QueryPerformanceCounter(&m_begin);
#else
    clock_gettime(CLOCK_REALTIME, &m_begin);
#endif
	}

	void end()
	{
#if _WIN32
		QueryPerformanceCounter(&m_end);
		QueryPerformanceFrequency(&m_freq);
#else
    clock_gettime(CLOCK_REALTIME, &m_end);
#endif
	}

	double getMilliseconds() const
	{
#if _WIN32
		return double((m_end.QuadPart - m_begin.QuadPart) * 1000) / double(m_freq.QuadPart);
#else
    return double((m_end.tv_nsec - m_end.tv_nsec)/1e9);
#endif
	}

private:
#if _WIN32
	LARGE_INTEGER m_begin;
	LARGE_INTEGER m_end;
	LARGE_INTEGER m_freq;
#else
  timespec m_begin;
  timespec m_end;
#endif
};

void PrintCreationTimes()
{
	uint32_t count = 1000;

	std::cout << "UVSphere" << std::endl;
	for (uint32_t i = 0; i < 22; ++i)
	{
		Timer timer;
		timer.start();
		Mesh m;
		for (uint32_t j = 0; j < count; ++j)
		{
			m.clear();
			UVSphere(2 * (i + 2), i + 2, m);
		}
		timer.end();
		std::cout << m.triangleCount() << " (x" << count << "): " << timer.getMilliseconds() << "ms" << std::endl;
	}

	std::cout << "Normalized Cube" << std::endl;
	for (uint32_t i = 1; i < 14; ++i)
	{
		Timer timer;
		timer.start();
		Mesh m;
		for (uint32_t j = 0; j < count; ++j)
		{
			m.clear();
			NormalizedCube(i, m);
		}
		timer.end();
		std::cout << m.triangleCount() << " (x" << count << "): " << timer.getMilliseconds() << "ms" << std::endl;
	}

	std::cout << "Spherified Cube" << std::endl;
	for (uint32_t i = 1; i < 14; ++i)
	{
		Timer timer;
		timer.start();
		Mesh m;
		for (uint32_t j = 0; j < count; ++j)
		{
			m.clear();
			SpherifiedCube(i, m);
		}
		timer.end();
		std::cout << m.triangleCount() << " (x" << count << "): " << timer.getMilliseconds() << "ms" << std::endl;
	}

	std::cout << "Icosahedron" << std::endl;
	{
		Mesh meshes[5];
		for (uint32_t i = 0; i < 5; ++i)
		{
			Timer timer;
			timer.start();
			for (uint32_t j = 0; j < count; ++j)
			{
				for (auto &m : meshes) m.clear();
				Icosahedron(meshes[0]);
				for (uint32_t k = 0; k < i; ++k) SubdivideMesh(meshes[k], meshes[k+1]);
			}
			timer.end();
			std::cout << meshes[i].triangleCount() << " (x" << count << "): " << timer.getMilliseconds() << "ms" << std::endl;
		}
	}
}

int main(int argc, char *argv[])
{
#if EXPORT_OBJS
	ExportMeshes();
#endif

#if EXPORT_IMAGES
	ExportErrorImages();
#endif

#if PRINT_STATS
	PrintStats();
#endif

#if PRINT_CREATION_TIMES
	PrintCreationTimes();
#endif

	return 0;
}
