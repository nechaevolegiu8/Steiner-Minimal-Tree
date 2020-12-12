#ifndef H_DELAUNAY
#define H_DELAUNAY

#include "vector2.h"
#include "edge.h"
#include "triangle.h"

#include <vector>
#include <algorithm>
#include <omp.h>

template <class T>
class Delaunay
{

public:
	using TriangleType = Triangle<T>;
	using EdgeType = Edge<T>;
	using VertexType = Vector2<T>;

	const std::vector<TriangleType>& Load(char* Path)
	{
		std::vector<Vector2<float>> points;
		std::ifstream file(Path, std::ios_base::in);
		if (!file.is_open()) throw "Cant open file / Wrong way";

		int size = 0;
		while (!file.eof())
		{
			file.getline(new char[255], 255);
			size++;
		}
		file.seekg(0);

		for (int i = 0; i < size; i++)
		{
			int a, b; file >> a; file >> b;
			points.push_back(Vector2<float>(a, b));
		}
		file.seekg(0);

		return triangulate(points);
	}

	const std::vector<TriangleType>& triangulate(std::vector<VertexType> &vertices)
	{	
		
		// Store the vertices localy
		_vertices = vertices;

		// Determinate the super triangle
		float minX = vertices[0].x;
		float minY = vertices[0].y;
		float maxX = minX;
		float maxY = minY;

		for (std::size_t i = 0; i < vertices.size(); ++i)
		{
			if (vertices[i].x < minX) minX = vertices[i].x;
			if (vertices[i].y < minY) minY = vertices[i].y;
			if (vertices[i].x > maxX) maxX = vertices[i].x;
			if (vertices[i].y > maxY) maxY = vertices[i].y;
		}

		float dx = maxX - minX;
		float dy = maxY - minY;
		float deltaMax = std::max(dx, dy);
		float midx = (minX + maxX) / 2.f;
		float midy = (minY + maxY) / 2.f;

		VertexType p1(midx - 20 * deltaMax, midy - deltaMax);
		VertexType p2(midx, midy + 20 * deltaMax);
		VertexType p3(midx + 20 * deltaMax, midy - deltaMax);

		// Create a list of triangles, and add the supertriangle in it
		_triangles.push_back(TriangleType(p1, p2, p3));	

		for (auto p = begin(vertices); p != end(vertices); p++)
		{
			std::vector<EdgeType> polygon;
			
			for (auto t = begin(_triangles); t != end(_triangles); t++)
			{
				// Processing
				if (t->circumCircleContains(*p))
				{
					// Pushing bad triangle
					t->isBad = true;
					polygon.push_back(t->e1);
					polygon.push_back(t->e2);
					polygon.push_back(t->e3);
				}
				else
				{
					//std::cout << " does not contains " << *p << " in his circum center" << std::endl;
				}
			}

			_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [](TriangleType &t) {
				return t.isBad;
			}), end(_triangles));

			for (auto e1 = begin(polygon); e1 != end(polygon); e1++)
			{
				for (auto e2 = begin(polygon); e2 != end(polygon); e2++)
				{
					if (e1 == e2)
						continue;

					if (*e1 == *e2)
					{
						e1->isBad = true;
						e2->isBad = true;
					}
				}
			}

			polygon.erase(std::remove_if(begin(polygon), end(polygon), [](EdgeType &e) {
				return e.isBad;
			}), end(polygon));

			for (auto e = begin(polygon); e != end(polygon); e++)
				_triangles.push_back(TriangleType(e->p1, e->p2, *p));

		}

		_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [p1, p2, p3](TriangleType &t) {
			return t.containsVertex(p1) || t.containsVertex(p2) || t.containsVertex(p3);
		}), end(_triangles));

		for (auto t = begin(_triangles); t != end(_triangles); t++)
		{
			_edges.push_back(t->e1);
			_edges.push_back(t->e2);
			_edges.push_back(t->e3);
		}

		return _triangles;
	}

	const std::vector<TriangleType>& getTriangles() const { return _triangles; };
	const std::vector<EdgeType>& getEdges() const { return _edges; };
	const std::vector<VertexType>& getVertices() const { return _vertices; };

private:
	std::vector<TriangleType> _triangles;
	std::vector<EdgeType> _edges;
	std::vector<VertexType> _vertices;
};

#endif
