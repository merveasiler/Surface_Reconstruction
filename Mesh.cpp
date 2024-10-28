#include "Mesh.h"
//#include "KdTree.h"
#include "MCubes.h"

void Mesh::calculateVertexNormal(int id) {
	Vertex* v = verts[id];

	// take weighted average of normals
	float* normal = new float[3];
	normal[0] = 0; normal[1] = 0; normal[2] = 0;
	for (int i=0; i < v->triList.size(); i++) {
		Triangle* t = tris[ v->triList[i] ];
		float weight = 1;// 1 / t->area;
		normal[0] += t->normal[0] * weight;
		normal[1] += t->normal[1] * weight;
		normal[2] += t->normal[2] * weight;
	}

	// normalize
	float length = v->triList.size();
	normal[0] = normal[0]/length; normal[1] = normal[1]/length; normal[2] = normal[2]/length;
	length = sqrt(pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2));
	if (length > 0)
		normal[0] = normal[0]/length; normal[1] = normal[1]/length; normal[2] = normal[2]/length;
	v->normal = normal;

}

void Mesh::constructGrid(float densityNoise) {

	float gridEdge = maxZ-minZ;
	if (maxY-minY > gridEdge)
		gridEdge = maxY-minY;
	if (maxX-minX > gridEdge)
		gridEdge = maxX-minX;
	float cellEdgeLength = gridEdge / 50;
	int numOfCells = gridEdge/cellEdgeLength;

	for (int i=0; i<verts.size(); i++)
		calculateVertexNormal(i);

	// construct & classify grid vertices
	XYZ grid[numOfCells+3][numOfCells+3][numOfCells+3];
	for (int i=0; i<numOfCells+3; i++) {
		float xCoord = minX + (i-1)*cellEdgeLength;
		for (int j=0; j<numOfCells+3; j++) {
			float yCoord = minY +(j-1)*cellEdgeLength;
			for (int k=0; k<numOfCells+3; k++) {
				float zCoord = minZ + (k-1)*cellEdgeLength;

				// define grid vertex coordinates
				(grid[i][j][k]).x = xCoord;
				(grid[i][j][k]).y = yCoord;
				(grid[i][j][k]).z = zCoord;

				// find the nearest mesh vertex
				float min_dist;
				int nearest_id;
				for (int s=0; s<verts.size(); s++) {
					Vertex* vx = verts[s];
					float dist = sqrt( pow(vx->coords[0]-xCoord, 2) + 
							   pow(vx->coords[1]-yCoord, 2) + 
							   pow(vx->coords[2]-zCoord, 2) );
					if (s==0) { min_dist = dist; nearest_id = s; continue; }
					if (dist < min_dist) { min_dist = dist; nearest_id = s; }
				}

				// calculate signed distance function (sdf)
				if (min_dist < densityNoise*cellEdgeLength) {
					Vertex* v = verts[nearest_id];
					float* vertNormal = v->normal;
					float distVect[3];
					distVect[0] = (xCoord-(v->coords[0])); 
					distVect[1] = (yCoord-(v->coords[1])); 
					distVect[2] = (zCoord-(v->coords[2]));
					float distLength = sqrt(pow(distVect[0],2)+pow(distVect[1],2)+pow(distVect[2],2));
					float sdf = (distVect[0]/distLength) * vertNormal[0] + 
						    (distVect[1]/distLength) * vertNormal[1] +
						    (distVect[2]/distLength) * vertNormal[2];
					(grid[i][j][k]).sdf = sdf;
				}
				else	(grid[i][j][k]).sdf = 1;
			}
		}
	}

	// execute Marching Cubes Algorithm
	for (int i=0; i<numOfCells+2; i++) {
		for (int j=0; j<numOfCells+2; j++) {
			for (int k=0; k<numOfCells+2; k++) {
				GRIDCELL cell;
				cell.p[0] = grid[i][j][k];
				cell.val[0] = (grid[i][j][k]).sdf;
				cell.p[1] = grid[i+1][j][k];
				cell.val[1] = (grid[i+1][j][k]).sdf;
				cell.p[2] = grid[i+1][j][k+1];
				cell.val[2] = (grid[i+1][j][k+1]).sdf;
				cell.p[3] = grid[i][j][k+1];
				cell.val[3] = (grid[i][j][k+1]).sdf;
				cell.p[4] = grid[i][j+1][k];
				cell.val[4] = (grid[i][j+1][k]).sdf;
				cell.p[5] = grid[i+1][j+1][k];
				cell.val[5] = (grid[i+1][j+1][k]).sdf;
				cell.p[6] = grid[i+1][j+1][k+1];
				cell.val[6] = (grid[i+1][j+1][k+1]).sdf;
				cell.p[7] = grid[i][j+1][k+1];
				cell.val[7] = (grid[i][j+1][k+1]).sdf;

				TRIANGLE triangles[5];
				int number = Polygonise(cell, 0, triangles);
				if (number>0) {
					for (int t=0; t<number; t++) { // triangles
						int n = estimatedVerts.size();
						for (int v=0; v<3; v++) { // triangle vertices
							float* coords = new float[3];
							coords[0] = triangles[t].p[v].x;
							coords[1] = triangles[t].p[v].y;
							coords[2] = triangles[t].p[v].z;
							Vertex* vertex = new Vertex(n+v, coords);
							estimatedVerts.push_back(vertex);
						}
						Triangle* triangle = new Triangle(estimatedTris.size(), n, n+1, n+2);
						estimatedTris.push_back(triangle);
					}
				}
			}
		}
	}
}

void Mesh::loadOff(const char* name)
{
	int nVerts, nTris, n, id, i = 0;
	float x, y, z;
	string line;

	ifstream myFile (name);
	if (myFile.is_open())
	{
		getline (myFile,line);	// read "OFF"
		myFile >> nVerts >> nTris >> n;

		while (i++ < nVerts)
		{
			myFile >> x >> y >> z;
			addVertex(x, y, z);
		}

		i=0;
		while (i++ <nTris)
		{
			myFile >> id >> x >> y >> z;
			addTriangle((int) x, (int) y, (int) z);
		}
		myFile.close();
	}
}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = verts.size();
	float* c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;

	verts.push_back( new Vertex(idx, c) );

	// find max/min coordinates
	if (idx == 0) {
		maxX = x; minX = x;
		maxY = y; minY = y;
		maxZ = z; minZ = z;
	}
	else {
		if (x>maxX) maxX = x;
		else if (x<minX) minX = x;
		else;
		
		if (y>maxY) maxY = y;
		else if (y<minY) minY = y;
		else;

		if (z>maxZ) maxZ = z;
		else if (z<minZ) minZ = z;
		else;
	}
}


void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	Triangle* t = new Triangle(idx, v1, v2, v3);
	t->computeCentNormArea(verts[v1], verts[v2], verts[v3]);
	tris.push_back(t);

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	if (! makeVertsNeighbor(v1, v2, idx) )
		addEdge(v1, v2, idx);

	if (! makeVertsNeighbor(v1, v3, idx) )
		addEdge(v1, v3, idx);

	if (! makeVertsNeighbor(v2, v3, idx) )
		addEdge(v2, v3, idx);

}

bool Mesh::makeVertsNeighbor(int v1i, int v2i, int t_id)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < verts[v1i]->edgeList.size(); i++) {
		Edge* edge = edges[verts[v1i]->edgeList[i]];
		if (edge->v1i == v2i || edge->v2i == v2i) {

			Triangle* tNew = tris[t_id];
			tNew->edgeList.push_back(edge->idx);

			for (int j=0; j < edge->triList.size(); j++) {
				Triangle* tOld = tris[edge->triList[j]];
				tOld->triList.push_back(tNew->idx);
				tNew->triList.push_back(tOld->idx);
			}

			edge->triList.push_back(t_id);
			return true;
		}
	}


	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return false;
}

void Mesh::addEdge(int v1, int v2, int t_id)
{
	int idx = edges.size();
	edges.push_back( new Edge(idx, v1, v2) );
	tris[t_id]->edgeList.push_back(idx);
	edges[idx]->triList.push_back(t_id);
	edges[idx]->computeLength(verts[v1], verts[v2]);

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

Mesh::~Mesh() {
	int i;
	for (i=0; i<tris.size(); i++) {
		delete[] tris[i]->center;
		delete[] tris[i]->normal;
		(tris[i]->triList).clear();
		(tris[i]->edgeList).clear();
		delete tris[i];
	}

	for (i=0; i<edges.size(); i++) {
		(edges[i]->triList).clear();
		delete edges[i];
	}

	for (i=0; i<verts.size(); i++) {
		(verts[i]->vertList).clear();
		(verts[i]->triList).clear();
		(verts[i]->edgeList).clear();
		delete[] verts[i]->coords;
		delete[] verts[i]->normal;
		delete verts[i];
	}

	for (i=0; i<estimatedVerts.size(); i++) {
		delete[] estimatedVerts[i]->coords;
		delete estimatedVerts[i];
	}

	for (i=0; i<estimatedTris.size(); i++)
		delete estimatedTris[i];

	verts.clear();
	edges.clear();
	tris.clear();
	estimatedVerts.clear();
	estimatedTris.clear();
}

float computeDet(float* column1, float* column2, float* column3) {
	float part1 = column1[0] * (column2[1]*column3[2] - column2[2]*column3[1]);
	float part2 = -column2[0] * (column1[1]*column3[2] - column1[2]*column3[1]);
	float part3 = column3[0] * (column1[1]*column2[2] - column1[2]*column2[1]);
	return part1+part2+part3;
}

float* Mesh::intersectRayTriangle(float* start, float* dir) {
	
	float* intersection = new float[5]; // the intersection point [0,1,2], triangle id [3] and tVal [4], resp.
	intersection[0] = 0, intersection[1] = 0, intersection[2] = 0, intersection[3] = -1, intersection[4] = -1;

	for (unsigned int i=0; i<tris.size(); i++) {

		Triangle* t = tris[i];
		float* a = verts[t->getVi(0)]->coords;
		float* b = verts[t->getVi(1)]->coords;
		float* c = verts[t->getVi(2)]->coords;

		float resultColumn[3] = {a[0]-start[0], a[1]-start[1], a[2]-start[2]};
		float column1[3] = {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
		float column2[3] = {a[0]-c[0], a[1]-c[1], a[2]-c[2]};

		float detA = computeDet(column1, column2, dir);
		float detBeta = computeDet(resultColumn, column2, dir);
		float detGama = computeDet(column1, resultColumn, dir);
		float det_t = computeDet(column1, column2, resultColumn);

		float beta = detBeta / detA;
		float gama = detGama / detA;
		float alpha = 1 - (beta+gama);
		float tVal = det_t / detA;
		float epsilon = 0.000001;

		if (alpha >= epsilon && alpha <= 1 && beta >= epsilon && beta <= 1 && gama >= epsilon && gama <= 1 && tVal>0) {
			intersection[0] = start[0] + dir[0]*tVal;
			intersection[1] = start[1] + dir[1]*tVal;
			intersection[2] = start[2] + dir[2]*tVal;
			intersection[3] = i;
			intersection[4] = tVal;
			return intersection;
		}
	}

	return intersection;
}

void Triangle::computeCentNormArea(Vertex* v1, Vertex* v2, Vertex* v3) {
	// compute center
	float x = (v1->coords[0] + v2->coords[0])/2;
	float y = (v1->coords[1] + v2->coords[1])/2;
	float z = (v1->coords[2] + v2->coords[2])/2;

	center = new float[3];
	center[0] = v3->coords[0] + (x - v3->coords[0])*2/3.0;
	center[1] = v3->coords[1] + (y - v3->coords[1])*2/3.0;
	center[2] = v3->coords[2] + (z - v3->coords[2])*2/3.0;

	// compute normal
	float x1 = v2->coords[0] - v1->coords[0];
	float y1 = v2->coords[1] - v1->coords[1];
	float z1 = v2->coords[2] - v1->coords[2];
	float x2 = v3->coords[0] - v1->coords[0];
	float y2 = v3->coords[1] - v1->coords[1];
	float z2 = v3->coords[2] - v1->coords[2];

	x = y1*z2 - y2*z1;
	y = z1*x2 - z2*x1;
	z = x1*y2 - x2*y1;
	float length = sqrt(x*x + y*y + z*z);

	normal = new float[3];
	normal[0] = x/length;
	normal[1] = y/length;
	normal[2] = z/length;

	// compute area
	area = length / 2;
}

int Triangle::findCommonEdge(Triangle* neighbor) {
	for (int i=0; i < edgeList.size(); i++) {
		int edgeId = edgeList[i];
		for (int j=0; j < neighbor->edgeList.size(); j++)
			if (edgeId == neighbor->edgeList[j])
				return edgeId;
	}
}

