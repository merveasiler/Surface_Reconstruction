#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

using Eigen::MatrixXd;
using namespace std;

struct gridVertex {
	float coords[3];
	float sdf;
};

struct Vertex
{
	float* coords, * normal; //3d coordinates etc
	int idx; //who am i; verts[idx]

	vector< int > vertList; //adj vertices;
	vector< int > triList; 
	vector< int > edgeList; 
	
	Vertex(int i, float* c) : idx(i), coords(c) {};
	float getCoord(int i) const {return coords[i];} ;
	float* getCoords() const {return coords;} ;
};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	vector< int > triList; //adj triangles (at max 2)
	float length;
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2) {};

	void computeLength(Vertex* v1, Vertex* v2)
	{
		float* coords1 = v1->coords;
		float* coords2 = v2->coords;
		length = sqrt(pow(coords1[0]-coords2[0], 2) + pow(coords1[1]-coords2[1], 2) + pow(coords1[2]-coords2[2], 2));
	}
};

struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
	float* center, * normal;
	float area;

	vector< int > edgeList;
	vector< int > triList; //adj triangles (at max 3)

	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {} ;
	int getVi(int i) const {if (i==0) return v1i; else if (i==1) return v2i; else return v3i;} ;
	void computeCentNormArea(Vertex* v1, Vertex* v2, Vertex* v3);
	int findCommonEdge(Triangle* neighbor);
};

class Mesh
{
private:
	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;

	float minX, maxX, minY, maxY, minZ, maxZ;
	vector< Triangle* > estimatedTris;
	vector< Vertex* > estimatedVerts;

	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2, int t_id);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i, int t_id);

public:

	Mesh() {} ;
	~Mesh();
	const Vertex* getVertex(int id) {return verts[id];} ;
	const Triangle* getTriangle(int id) {return tris[id];} ;
	const Edge* getEdge(int id) {return edges[id];} ;
	int getNumOfVerts() {return verts.size();} ;
	int getNumOfEdges() {return edges.size();} ;
	int getNumOfTris() {return tris.size();} ;
	void loadOff(const char* name);
	float* intersectRayTriangle(float* start, float* dir);
	void calculateVertexNormal(int id);
	void constructGrid(float densityNoise);
	int getNumOfEstVerts() {return estimatedVerts.size();} ;
	int getNumOfEstTris() {return estimatedTris.size();} ;
	const Vertex* getEstVertex(int id) {return estimatedVerts[id];} ;
	const Triangle* getEstTriangle(int id) {return estimatedTris[id];} ;
};

