#ifndef H_STEINER
#define H_STEINER
#include <math.h>
#include "vector2.h"
#include "triangle.h"
#define PI 3.14159265 // Mysterious number

template <class T>
class Steiner
{
public:
	using TriangleType = Triangle<T>;
	using VertexType = Vector2<T>;

	const std::vector<VertexType>& additionalVertices(std::vector<TriangleType> &triangles)
	{
		std::vector<VertexType> additionalvertices;

		omp_set_num_threads(2);
		#pragma omp parallel for // Using OpenMP
		for (int e1 = 0; e1 < triangles.size(); e1++)
		{
			VertexType Vertex, RadiusVertex, SteinerVertex;
			int n = largestAngle(triangles[e1].p1, triangles[e1].p2, triangles[e1].p3);
			
			if (n == 1)
			{
				Vertex = findThirdVertex(triangles[e1].p2, triangles[e1].p1, triangles[e1].p3);
				RadiusVertex = findCenterVertex(triangles[e1].p2, triangles[e1].p1, triangles[e1].p3);
				SteinerVertex = findSteinerVertex(triangles[e1].p1, Vertex, RadiusVertex, triangles[e1].p2, triangles[e1].p3);
			}
			
			if (n == 2)
			{
				Vertex = findThirdVertex(triangles[e1].p1, triangles[e1].p2, triangles[e1].p3);
				RadiusVertex = findCenterVertex(triangles[e1].p1, triangles[e1].p2, triangles[e1].p3);
				SteinerVertex = findSteinerVertex(triangles[e1].p2, Vertex, RadiusVertex, triangles[e1].p1, triangles[e1].p3);
			}
			
			if (n == 3)
			{
				Vertex = findThirdVertex(triangles[e1].p1, triangles[e1].p3, triangles[e1].p2);
				RadiusVertex = findCenterVertex(triangles[e1].p1, triangles[e1].p3, triangles[e1].p2);
				SteinerVertex = findSteinerVertex(triangles[e1].p3, Vertex, RadiusVertex, triangles[e1].p1, triangles[e1].p2);
			}
			
			if (n == 4)
			{
				//	std::cout << "Triangle #" << e1 + 1 << std::endl << "No steiner vertex (Max angle > 120*) " << std::endl;
				std::cout << ".";
			}
			
			else 
			{
				// std::cout << "Triangle #" << e1 + 1 << std::endl << "Steiner vertex: " << SteinerVertex << std::endl;
				std::cout << ".";
				_vertices.push_back(SteinerVertex);
			}
		}
		
		return _vertices;
	}

	// Function to find distance between vertices
	const float l(VertexType &v1, VertexType &v2)
	{
		return sqrt(pow((v2.x - v1.x), 2) + pow((v2.y - v1.y), 2));
	}

	// Function to find which vertex has largest angle
	const int largestAngle(VertexType &v1, VertexType &v2, VertexType &v3)
	{
		float max = 0; int which;

		float A = acos((pow(l(v1, v2), 2) + pow(l(v1, v3), 2) - pow(l(v2, v3), 2))
			/ (2 * l(v1, v2) * l(v1, v3))) * 180.0 / PI;
		float B = acos((pow(l(v1, v2), 2) + pow(l(v2, v3), 2) - pow(l(v1, v3), 2))
			/ (2 * l(v1, v2) * l(v2, v3))) * 180.0 / PI;
		float C = acos((pow(l(v1, v3), 2) + pow(l(v2, v3), 2) - pow(l(v1, v2), 2))
			/ (2 * l(v1, v3) * l(v2, v3))) * 180.0 / PI;

		if (A > max) { max = A; which = 1; }
		if (B > max) { max = B; which = 2; }
		if (C > max) { max = C; which = 3; }
		if (max > 120) return 4;
		else return which;
	}

	// Function to find third vertex of new triangle with equal sides
	const VertexType& findThirdVertex(VertexType &p1, VertexType &p2, VertexType &p3)
	{
		float x2 = (4 * pow(p1.x, 3) - 4 * pow(p1.x, 2) * p3.x + sqrt(pow((-4 * pow(p1.x, 3) + 4 * pow(p1.x, 2) * p3.x + 4 * p1.x * pow(p3.x, 2)
			- 4 * p1.x * pow(p1.y, 2) + 8 * p1.x * p1.y * p3.y - 4 * p1.x * pow(p3.y, 2) - 4 * pow(p3.x, 3) - 4 * p3.x * pow(p1.y, 2) + 8 * p3.x
			* p1.y * p3.y - 4 * p3.x * pow(p3.y, 2)), 2) - 4 * (4 * pow(p1.x, 2) - 8 * p1.x * p3.x + 4 * pow(p3.x, 2) + 4 * pow(p1.y, 2) - 8 * p1.y
			* p3.y + 4 * pow(p3.y, 2)) * (pow(p1.x, 4) - 2 * pow(p1.x, 2) * pow(p3.x, 2) + 2 * pow(p1.x, 2) * pow(p1.y, 2) - 4 * pow(p1.x, 2) * p1.y
			* p3.y + 2 * pow(p1.x, 2) * pow(p3.y, 2) + pow(p3.x, 4) + 2 * pow(p3.x, 2) * pow(p1.y, 2) - 4 * pow(p3.x, 2) * p1.y * p3.y + 2 * pow(p3.x, 2)
			* pow(p3.y, 2) + pow(p1.y, 4) - 4 * pow(p1.y, 3) * p3.y + 6 * pow(p1.y, 2) * pow(p3.y, 2) - 4 * pow(p1.y, 2) * pow(l(p1, p3), 2) - 4 * p1.y
			* pow(p3.y, 3) + 8 * p1.y * p3.y * pow(l(p1, p3), 2) + pow(p3.y, 4) - 4 * pow(p3.y, 2) * pow(l(p1, p3), 2))) - 4 * p1.x * pow(p3.x, 2) + 4
			* p1.x * pow(p1.y, 2) - 8 * p1.x * p1.y * p3.y + 4 * p1.x * pow(p3.y, 2) + 4 * pow(p3.x, 3) + 4 * p3.x * pow(p1.y, 2) - 8 * p3.x * p1.y
			* p3.y + 4 * p3.x * pow(p3.y, 2)) / (2 * (4 * pow(p1.x, 2) - 8 * p1.x * p3.x + 4 * pow(p3.x, 2) + 4 * pow(p1.y, 2) - 8 * p1.y * p3.y
			+ 4 * pow(p3.y, 2)));
		
		float y2 = p1.y - sqrt(-pow(p1.x, 2) + 2 * p1.x * x2 - pow(x2, 2) + pow(l(p1, p3), 2));
		VertexType point1(x2, y2);
		
		if (isEqual(l(point1, p1), l(point1, p3)) != 1) {
			y2 = sqrt(-pow(p1.x, 2) + 2 * p1.x * x2 - pow(x2, 2) + pow(l(p1, p3), 2)) + p1.y;
			VertexType point2(x2, y2); point1 = point2;
		}

		x2 = (4 * pow(p1.x, 3) - 4 * pow(p1.x, 2) * p3.x - sqrt(pow((-4 * pow(p1.x, 3) + 4 * pow(p1.x, 2) * p3.x + 4 * p1.x * pow(p3.x, 2)
			- 4 * p1.x * pow(p1.y, 2) + 8 * p1.x * p1.y * p3.y - 4 * p1.x * pow(p3.y, 2) - 4 * pow(p3.x, 3) - 4 * p3.x * pow(p1.y, 2) + 8 * p3.x
			* p1.y * p3.y - 4 * p3.x * pow(p3.y, 2)), 2) - 4 * (4 * pow(p1.x, 2) - 8 * p1.x * p3.x + 4 * pow(p3.x, 2) + 4 * pow(p1.y, 2) - 8 * p1.y
			* p3.y + 4 * pow(p3.y, 2)) * (pow(p1.x, 4) - 2 * pow(p1.x, 2) * pow(p3.x, 2) + 2 * pow(p1.x, 2) * pow(p1.y, 2) - 4 * pow(p1.x, 2) * p1.y
			* p3.y + 2 * pow(p1.x, 2) * pow(p3.y, 2) + pow(p3.x, 4) + 2 * pow(p3.x, 2) * pow(p1.y, 2) - 4 * pow(p3.x, 2) * p1.y * p3.y + 2 * pow(p3.x, 2)
			* pow(p3.y, 2) + pow(p1.y, 4) - 4 * pow(p1.y, 3) * p3.y + 6 * pow(p1.y, 2) * pow(p3.y, 2) - 4 * pow(p1.y, 2) * pow(l(p1, p3), 2) - 4 * p1.y
			* pow(p3.y, 3) + 8 * p1.y * p3.y * pow(l(p1, p3), 2) + pow(p3.y, 4) - 4 * pow(p3.y, 2) * pow(l(p1, p3), 2))) - 4 * p1.x * pow(p3.x, 2) + 4
			* p1.x * pow(p1.y, 2) - 8 * p1.x * p1.y * p3.y + 4 * p1.x * pow(p3.y, 2) + 4 * pow(p3.x, 3) + 4 * p3.x * pow(p1.y, 2) - 8 * p3.x * p1.y
			* p3.y + 4 * p3.x * pow(p3.y, 2)) / (2 * (4 * pow(p1.x, 2) - 8 * p1.x * p3.x + 4 * pow(p3.x, 2) + 4 * pow(p1.y, 2) - 8 * p1.y * p3.y
			+ 4 * pow(p3.y, 2)));

		y2 = sqrt(-pow(p1.x, 2) + 2 * p1.x * x2 - pow(x2, 2) + pow(l(p1, p3), 2)) + p1.y;
		VertexType point3(x2, y2);
		
		if (isEqual(l(point3, p1), l(point3, p3)) != 1) {
			y2 = p1.y - sqrt(-pow(p1.x, 2) + 2 * p1.x * x2 - pow(x2, 2) + pow(l(p1, p3), 2));
			VertexType point4(x2, y2); point3 = point4;
		}

		if (l(p2, point1) > l(p2, point3)) return point1;
		else return point3;
	}

	// Function to find center of new triangle with equal sides
	const VertexType& findCenterVertex(VertexType &p1, VertexType &p2, VertexType &p3)
	{
		float x2 = (12 * pow(p1.x, 3) - 12 * pow(p1.x, 2) * p3.x - sqrt(pow((-12 * pow(p1.x, 3) + 12 * pow(p1.x, 2) * p3.x + 12 * p1.x * pow(p3.x, 2)
			- 12 * p1.x * pow(p1.y, 2) + 24 * p1.x * p1.y * p3.y - 12 * p1.x * pow(p3.y, 2) - 12 * pow(p3.x, 3) - 12 * p3.x * pow(p1.y, 2) + 24 * p3.x
			* p1.y * p3.y - 12 * p3.x * pow(p3.y, 2)), 2) - 4 * (12 * pow(p1.x, 2) - 24 * p1.x * p3.x + 12 * pow(p3.x, 2) + 12 * pow(p1.y, 2) - 24 * p1.y
			* p3.y + 12 * pow(p3.y, 2)) * (3 * pow(p1.x, 4) - 6 * pow(p1.x, 2) * pow(p3.x, 2) + 6 * pow(p1.x, 2) * pow(p1.y, 2) - 12 * pow(p1.x, 2) * p1.y
			* p3.y + 6 * pow(p1.x, 2) * pow(p3.y, 2) + 3 * pow(p3.x, 4) + 6 * pow(p3.x, 2) * pow(p1.y, 2) - 12 * pow(p3.x, 2) * p1.y * p3.y + 6 * pow(p3.x, 2)
			* pow(p3.y, 2) + 3 * pow(p1.y, 4) - 12 * pow(p1.y, 3) * p3.y + 18 * pow(p1.y, 2) * pow(p3.y, 2) - 4 * pow(p1.y, 2) * pow(l(p1, p3), 2)  - 12 * p1.y
			* pow(p3.y, 3) + 8 * p1.y * p3.y * pow(l(p1, p3), 2)  + 3 * pow(p3.y, 4) - 4 * pow(p3.y, 2) * pow(l(p1, p3), 2) )) - 12 * p1.x * pow(p3.x, 2) + 12
			* p1.x * pow(p1.y, 2) - 24 * p1.x * p1.y * p3.y + 12 * p1.x * pow(p3.y, 2) + 12 * pow(p3.x, 3) + 12 * p3.x * pow(p1.y, 2) - 24 * p3.x * p1.y
			* p3.y + 12 * p3.x * pow(p3.y, 2)) / (2 * (12 * pow(p1.x, 2) - 24 * p1.x * p3.x + 12 * pow(p3.x, 2) + 12 * pow(p1.y, 2) - 24 * p1.y * p3.y
			+ 12 * pow(p3.y, 2)));
		
		float y2 = sqrt(-3 * pow(p1.x, 2) + 6 * p1.x * x2 - 3 * pow(x2, 2) + pow(l(p1, p3), 2) ) / sqrt(3) + p1.y;
		VertexType point1(x2, y2);

		if (isEqual(l(point1, p1), l(point1, p3)) != 1) {
			y2 = p1.y - sqrt(-3 * pow(p1.x, 2) + 6 * p1.x * x2 - 3 * pow(x2, 2) + pow(l(p1, p3), 2)) / sqrt(3);
			VertexType point2(x2, y2); point1 = point2;
		}
		
		x2 = (12 * pow(p1.x, 3) - 12 * pow(p1.x, 2) * p3.x + sqrt(pow((-12 * pow(p1.x, 3) + 12 * pow(p1.x, 2) * p3.x + 12 * p1.x * pow(p3.x, 2)
			- 12 * p1.x * pow(p1.y, 2) + 24 * p1.x * p1.y * p3.y - 12 * p1.x * pow(p3.y, 2) - 12 * pow(p3.x, 3) - 12 * p3.x * pow(p1.y, 2) + 24
			* p3.x * p1.y * p3.y - 12 * p3.x * pow(p3.y, 2)), 2) - 4 * (12 * pow(p1.x, 2) - 24 * p1.x * p3.x + 12 * pow(p3.x, 2) + 12 * pow(p1.y, 2)
			- 24 * p1.y * p3.y + 12 * pow(p3.y, 2)) * (3 * pow(p1.x, 4) - 6 * pow(p1.x, 2) * pow(p3.x, 2) + 6 * pow(p1.x, 2) * pow(p1.y, 2) - 12
			* pow(p1.x, 2) * p1.y * p3.y + 6 * pow(p1.x, 2) * pow(p3.y, 2) + 3 * pow(p3.x, 4) + 6 * pow(p3.x, 2) * pow(p1.y, 2) - 12 * pow(p3.x, 2)
			* p1.y * p3.y + 6 * pow(p3.x, 2) * pow(p3.y, 2) + 3 * pow(p1.y, 4) - 12 * pow(p1.y, 3) * p3.y + 18 * pow(p1.y, 2) * pow(p3.y, 2) - 4
			* pow(p1.y, 2) * pow(l(p1, p3), 2)  - 12 * p1.y * pow(p3.y, 3) + 8 * p1.y * p3.y * pow(l(p1, p3), 2)  + 3 * pow(p3.y, 4) - 4 
			* pow(p3.y, 2) * pow(l(p1, p3), 2) )) - 12 * p1.x * pow(p3.x, 2) + 12 * p1.x * pow(p1.y, 2) - 24 * p1.x * p1.y * p3.y + 12
			* p1.x * pow(p3.y, 2) + 12 * pow(p3.x, 3) + 12 * p3.x * pow(p1.y, 2) - 24 * p3.x * p1.y * p3.y + 12 * p3.x * pow(p3.y, 2))
			/ (2 * (12 * pow(p1.x, 2) - 24 * p1.x * p3.x + 12 * pow(p3.x, 2) + 12 * pow(p1.y, 2) - 24 * p1.y * p3.y + 12 * pow(p3.y, 2)));
		
		y2 = p1.y - sqrt(-3 * pow(p1.x, 2) + 6 * p1.x * x2 - 3 * pow(x2, 2) + pow(l(p1, p3), 2) ) / sqrt(3);
		VertexType point3(x2, y2);
		
		if (isEqual(l(point3, p1), l(point3, p3)) != 1) {
			y2 = sqrt(-3 * pow(p1.x, 2) + 6 * p1.x * x2 - 3 * pow(x2, 2) + pow(l(p1, p3), 2)) / sqrt(3) + p1.y;
			VertexType point4(x2, y2); point3 = point4;
		}

		if (l(p2, point1) > l(p2, point3)) return point1;
		else return point3;
	}

	// Function to find steiner vertex in original triangle
	const VertexType& findSteinerVertex(VertexType &p2, VertexType &p3, VertexType &p4, VertexType &a, VertexType &c)
	{
		float y1 = (-sqrt(pow((-6 * pow(p2.x, 2) * p3.y + 6 * p2.x * p3.x * p2.y + 6 * p2.x * p3.x * p3.y - 6 * p2.x * p4.x
			* p2.y + 6 * p2.x * p4.x * p3.y - 6 * pow(p3.x, 2) * p2.y + 6 * p3.x * p4.x * p2.y - 6 * p3.x * p4.x * p3.y - 6
			* pow(p2.y, 2) * p4.y + 12 * p2.y * p3.y * p4.y - 6 * pow(p3.y, 2) * p4.y), 2) - 4 * (3 * pow(p2.x, 2) - 6 * p2.x
			* p3.x + 3 * pow(p3.x, 2) + 3 * pow(p2.y, 2) - 6 * p2.y * p3.y + 3 * pow(p3.y, 2)) * (3 * pow(p2.x, 2) * pow(p3.y, 2)
			- 6 * p2.x * p3.x * p2.y * p3.y + 6 * p2.x * p4.x * p2.y * p3.y - 6 * p2.x * p4.x * pow(p3.y, 2) + 3 * pow(p3.x, 2)
			* pow(p2.y, 2) - 6 * p3.x * p4.x * pow(p2.y, 2) + 6 * p3.x * p4.x * p2.y * p3.y + 3 * pow(p4.x, 2) * pow(p2.y, 2)
			- 6 * pow(p4.x, 2) * p2.y * p3.y + 3 * pow(p4.x, 2) * pow(p3.y, 2) + 3 * pow(p2.y, 2) * pow(p4.y, 2) - pow(p2.y, 2)
			* pow(l(a, c), 2) - 6 * p2.y * p3.y * pow(p4.y, 2) + 2 * p2.y * p3.y * pow(l(a, c), 2) + 3 * pow(p3.y, 2) * pow(p4.y, 2)
			- pow(p3.y, 2) * pow(l(a, c), 2))) + 6 * pow(p2.x, 2) * p3.y - 6 * p2.x * p3.x * p2.y - 6 * p2.x * p3.x * p3.y + 6 * p2.x
			* p4.x * p2.y - 6 * p2.x * p4.x * p3.y + 6 * pow(p3.x, 2) * p2.y - 6 * p3.x * p4.x * p2.y + 6 * p3.x * p4.x * p3.y + 6 
			* pow(p2.y, 2) * p4.y - 12 * p2.y * p3.y * p4.y + 6 * pow(p3.y, 2) * p4.y) / (2 * (3 * pow(p2.x, 2) - 6 * p2.x * p3.x 
			+ 3 * pow(p3.x, 2) + 3 * pow(p2.y, 2) - 6 * p2.y * p3.y + 3 * pow(p3.y, 2)));

		float x1 = p4.x - sqrt(abs(-3 * pow(y1, 2) + 6 * y1 * p4.y - 3 * pow(p4.y, 2) + pow(l(a, c), 2))) / sqrt(3);
		VertexType point1(x1, y1);

		x1 = p4.x + sqrt(abs(-3 * pow(y1, 2) + 6 * y1 * p4.y - 3 * pow(p4.y, 2) + pow(l(a, c), 2))) / sqrt(3);
		VertexType point2(x1, y1);

		y1 = (sqrt(pow((-6 * pow(p2.x, 2) * p3.y + 6 * p2.x * p3.x * p2.y + 6 * p2.x * p3.x * p3.y - 6 * p2.x * p4.x * p2.y + 6 
			* p2.x * p4.x * p3.y - 6 * pow(p3.x, 2) * p2.y + 6 * p3.x * p4.x * p2.y - 6 * p3.x * p4.x * p3.y - 6 * pow(p2.y, 2) * p4.y 
			+ 12 * p2.y * p3.y * p4.y - 6 * pow(p3.y, 2) * p4.y), 2) - 4 * (3 * pow(p2.x, 2) - 6 * p2.x * p3.x + 3 * pow(p3.x, 2) + 3
			* pow(p2.y, 2) - 6 * p2.y * p3.y + 3 * pow(p3.y, 2)) * (3 * pow(p2.x, 2) * pow(p3.y, 2) - 6 * p2.x * p3.x * p2.y * p3.y + 6
			* p2.x * p4.x * p2.y * p3.y - 6 * p2.x * p4.x * pow(p3.y, 2) + 3 * pow(p3.x, 2) * pow(p2.y, 2) - 6 * p3.x * p4.x * pow(p2.y, 2)
			+ 6 * p3.x * p4.x * p2.y * p3.y + 3 * pow(p4.x, 2) * pow(p2.y, 2) - 6 * pow(p4.x, 2) * p2.y * p3.y + 3 * pow(p4.x, 2) * pow(p3.y, 2)
			+ 3 * pow(p2.y, 2) * pow(p4.y, 2) - pow(p2.y, 2) * pow(l(a, c), 2) - 6 * p2.y * p3.y * pow(p4.y, 2) + 2 * p2.y * p3.y * pow(l(a, c), 2)
			+ 3 * pow(p3.y, 2) * pow(p4.y, 2) - pow(p3.y, 2) * pow(l(a, c), 2))) + 6 * pow(p2.x, 2) * p3.y - 6 * p2.x * p3.x * p2.y - 6 * p2.x * p3.x
			* p3.y + 6 * p2.x * p4.x * p2.y - 6 * p2.x * p4.x * p3.y + 6 * pow(p3.x, 2) * p2.y - 6 * p3.x * p4.x * p2.y + 6 * p3.x * p4.x * p3.y + 6
			* pow(p2.y, 2) * p4.y - 12 * p2.y * p3.y * p4.y + 6 * pow(p3.y, 2) * p4.y) / (2 * (3 * pow(p2.x, 2) - 6 * p2.x * p3.x + 3 * pow(p3.x, 2)
			+ 3 * pow(p2.y, 2) - 6 * p2.y * p3.y + 3 * pow(p3.y, 2)));
		
		x1 = p4.x - sqrt(abs(-3 * pow(y1, 2) + 6 * y1 * p4.y - 3 * pow(p4.y, 2) + pow(l(a, c), 2))) / sqrt(3);
		VertexType point3(x1, y1);

		x1 = p4.x + sqrt(abs(-3 * pow(y1, 2) + 6 * y1 * p4.y - 3 * pow(p4.y, 2) + pow(l(a, c), 2))) / sqrt(3);
		VertexType point4(x1, y1);
		
		int n = whichPoint(point1, point2, point3, point4, p2, a, c);
		
		if (n == 1) return point1;
		if (n == 2) return point2;
		if (n == 3) return point3;
		if (n == 4) return point4;
	}

	// Function to compare floats
	bool isEqual(float x, float y)
	{
		return std::fabs(x - y) < 0.001;
	}

	// Function to determine which point is point that we need
	int whichPoint(VertexType &p1, VertexType &p2, VertexType &p3, VertexType &p4, VertexType &a, VertexType &b, VertexType &c)
	{
		float min = l(p1, a) + l(p1, b) + l(p1, c); int n = 1;

		if ((l(p2, a) + l(p2, b) + l(p2, c)) < min) { min = l(p2, a) + l(p2, b) + l(p2, c); n = 2; }
		
		if ((l(p3, a) + l(p3, b) + l(p3, c)) < min) { min = l(p3, a) + l(p3, b) + l(p3, c); n = 3; }
		
		if ((l(p4, a) + l(p4, b) + l(p4, c)) < min) { min = l(p4, a) + l(p4, b) + l(p4, c); n = 4; }

		return n;
	}

private:
	std::vector<VertexType> _vertices;
	std::vector<TriangleType> _triangles;
};


#endif
