#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <iterator>
#include <omp.h>
#include <chrono>
//-----------
#include "vector2.h"
#include "triangle.h"
#include "delaunay.h"
#include "steiner.h"
#include "prim.h"

using namespace std::chrono;

int main()
{
	char path[255] = "files/good.dat";
	
	// Divide our graph into triangles (Delaunay triangulation) -(1)-
	
	auto start = high_resolution_clock::now(); // Time count start

	Delaunay<float> triangulation;
	std::vector<Triangle<float>> triangles = triangulation.Load(path);
	
	auto stop = high_resolution_clock::now(); // Time count stop
	auto duration1 = duration_cast<milliseconds>(stop - start); // Time count

	// Create additional vertices -------------------------------(2)-
	
	start = high_resolution_clock::now(); 
	
	Steiner<float> steiner;
	std::vector<Vector2<float>> steinerpoints = steiner.additionalVertices(triangles);

	stop = high_resolution_clock::now(); 
	auto duration2 = duration_cast<milliseconds>(stop - start); 

	// Find shortest path ---------------------------------------(3)-

	start = high_resolution_clock::now(); 

	std::vector<Vector2<float>> points = triangulation.getVertices();
	Prim<float> prim;
	
	float result = prim.shortestPath(points, steinerpoints); // Provides final solution

	stop = high_resolution_clock::now(); 
	auto duration3 = duration_cast<milliseconds>(stop - start); 

	// Show execution time for every part -----------------------(4)-

	std::cout << std::endl << "Time: " << std::endl << "Delaunay: " << duration1.count() << std::endl;
	std::cout << "Steiner:  " << duration2.count() << std::endl;
	std::cout << "Prim:     " << duration3.count() << std::endl;

	return 0;
}
