#ifndef H_PRIM
#define H_PRIM

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "vector2.h"
#include "edge.h"
#include "triangle.h"
#include "delaunay.h"

template <class T>
class Prim
{
public:
	using TriangleType = Triangle<T>;
	using VertexType = Vector2<T>;
	
	// Only for testing execution time of Prim algorithm, doesn't give final solution of SMT problem
	const float testTime(std::vector<VertexType>& vertices, std::vector<VertexType>& steinerpoints) {
		float min = FLT_MAX; int n;
		std::vector<float> results;
		std::vector<std::vector<int>> bin;

		std::vector<Vector2<float>> randPoints;

		for (int j = 0; j < vertices.size(); j++)
			randPoints.push_back(vertices[j]);

		std::vector<std::vector<float>> adjMatrix = getAdjMatrix(randPoints);
		return primMST(adjMatrix, 0);
	}

	const float shortestPath(std::vector<VertexType> &vertices, std::vector<VertexType> &steinerpoints)
	{
		float min = FLT_MAX; int n;
		std::vector<float> results;
		std::vector<std::vector<int>> bin; 
		
		// Lexicographic binary order for N Ferma points
		for (int i = 0; i < pow(2, steinerpoints.size()); i++)
		{
			char p[1024]; 
			_itoa_s(binary(i), p, 10);
			std::vector<int> temp;
			for (int j = 0; j < steinerpoints.size(); j++)
			{
				char s = p[j];
				temp.push_back(atoi(&s));
			}
			bin.push_back(temp);
		}
		
		// Bruteforce
		for (int i = 0; i < pow(2, steinerpoints.size()); i++)
		{
			std::vector<Vector2<float>> randPoints;
			for (auto &v : vertices) randPoints.push_back(v);
			for (int j = 0; j < steinerpoints.size(); j++)
			{
				if (bin[i][j] == 1)
					randPoints.push_back(steinerpoints[j]);
			}
			std::vector<std::vector<float>> adjMatrix = getAdjMatrix(randPoints);
			results.push_back(primMST(adjMatrix, 0));
		}
		
		// Finding best result
		for (int i = 0; i < results.size(); i++)
		{
			if (results[i] < min) { min = results[i]; n = i; }
		}

		// Show and return best result
		//std::cout << std::endl << "Best additional points:" << std::endl;
		std::vector<Vector2<float>> randPoints;
		
		for (auto &v : vertices) randPoints.push_back(v);
		for (int j = 0; j < steinerpoints.size(); j++)
		{
			if (bin[n][j] == 1) {
				randPoints.push_back(steinerpoints[j]);
				//std::cout << "x " << steinerpoints[j].x << " y " << steinerpoints[j].y << std::endl;
			}
		}
		std::vector<std::vector<float>> adjMatrix = getAdjMatrix(randPoints);

		std::cout << "Points, included in SMT: " << std::endl;
		for (int j = 0; j < randPoints.size(); j++)
		{
			std::cout << "#" << j+1 << "| x: " << randPoints[j].x << " | y: " << randPoints[j].y << " |" << std::endl;
		}

		primMST(adjMatrix, 1);
		return results[n];	
	}

	const std::vector<std::vector<float>> getAdjMatrix(std::vector<VertexType> &vertices)
	{
		std::vector<std::vector<float>> data;
		
		for (int i = 0; i < vertices.size(); i++)
		{
			std::vector<float> temp;
			for (int j = 0; j < vertices.size(); j++)
			{
				if (i == j)
				{ 
					temp.push_back(0); 
				}
				else
				{
					float f = l(vertices[i], vertices[j]);
					temp.push_back(f);
				}
			}
			data.push_back(std::vector<float>(temp));
		}

		return data;
	}

	const float l(VertexType &v1, VertexType &v2)
	{
		return sqrt(pow((v2.x - v1.x), 2) + pow((v2.y - v1.y), 2));
	}

	// Function to find the vertex with minimum key value, from the set 
	// of vertices not yet included in MST
	int minKey(std::vector<float> &key, std::vector<bool> &mstSet, std::vector<std::vector<float>> &graph)
	{
		// Initialize min value
		float min = FLT_MAX;
		int min_index = INT_MIN;

		omp_set_num_threads(2);
		#pragma omp parallel for // Using OpenMP
		for (int v = 0; v < graph.size(); v++) {
			if (mstSet[v] == false && key[v] < min)
				min = key[v], min_index = v;
		}

		return min_index;
	}

	// Function to print the constructed MST stored in parent[]
	const float Solution(std::vector<int> &parent, std::vector<std::vector<float>> &graph, int n) 
	{
		float summary = 0;
		if (n == 1)
		{
			printf("\nSolution: \n");
			printf("Path         Length\n");
			for (int i = 1; i < graph.size(); i++)
			{
				printf("#%d <-> #%d    %.2lf \n", parent[i] + 1, i + 1, graph[i][parent[i]]);
				summary = summary + graph[i][parent[i]];
			}
			printf("Summary: %.2lf \n", summary);
			return summary;
		}
		if (n == 0)
		{
			for (int i = 1; i < graph.size(); i++)
			{
				summary = summary + graph[i][parent[i]];
			}
			return summary;
		}
	}

	// Function to construct MST
	float primMST(std::vector<std::vector<float>> &graph, int n)
	{
		std::vector<int> parent; // Array to store constructed MST
		std::vector<float> key; // Key values used to pick minimum weight edge in cut
		std::vector<bool> mstSet; // To represent set of vertices not yet included in MST

		// Initialize all keys as INFINITE
		for (int i = 0; i < graph[0].size(); i++)
		{
			key.push_back(FLT_MAX);
			mstSet.push_back(false);
			parent.push_back(0);
		}

		// Always include first 1st vertex in MST
		key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
		parent[0] = -1; // First node is always root of MST

		for (int count = 0; count < graph.size() - 1; count++)
		{
			int u = 0;
			// Pick the minimum key vertex from the set of vertices
			// not yet included in MST 	
			u = minKey(key, mstSet, graph);
				// Add the picked vertex to the MST Set
				mstSet[u] = true;

				// Update the key only if graph[u][v] is smaller than key[v]
				for (int v = 0; v < graph.size(); v++)
					if (graph[u][v] && mstSet[v] == false && graph[u][v] < key[v])
						parent[v] = u, key[v] = graph[u][v];
		}

		// Print the constructed MST
		if (n == 1) return Solution(parent, graph, 1);
		if (n == 0) return Solution(parent, graph, 0);
	}
	
	int binary(int num)
	{
		int t = 0, d = 1;
		while (num)
		{
			t += (num % 2)*d;
			num = num / 2;
			d = d * 10;
		}
		return t;
	}

private:
	std::vector<VertexType> _vertices;
};

#endif 
