#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoCube.h>

#include "Mesh.h"
#include "Painter.h"

int main(int argc, char ** argv)
{
	QWidget * mainwin = SoQt::init(argc, argv, argv[0]);
	SoQtExaminerViewer * viewer = new SoQtExaminerViewer(mainwin);
	SoSeparator * root = new SoSeparator;
	root->ref();

	// get file name
	cout << "Enter the input file name: ";
	string inputFile;
	getline(cin, inputFile);

	// process file
	Painter* painter = new Painter();
	SoSeparator * res;
	Mesh* mesh = new Mesh();
	mesh->loadOff(inputFile.c_str());

	// get the operation type
	cout << "Enter the <density> + <noise> value: ";
	float densityNoise;
	cin >> densityNoise;

	// construct grid and do the 3D-reconstruction
	clock_t timeMeasure;
	timeMeasure = clock();
	mesh->constructGrid(densityNoise);
	timeMeasure = clock() - timeMeasure;
	cout << "It took " << ((float)timeMeasure)/CLOCKS_PER_SEC << " seconds to do the whole 3d-reconstruction operations.\n";

	res = painter->getShapeSep(mesh, false, false);
	//painter->drawTriangulation(res, mesh);

	// display on the screen
	root->addChild(res);
	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);	
	viewer->show();
	SoQt::show(mainwin);
	SoQt::mainLoop();
	delete viewer;
	root->unref();

	delete painter;
	delete mesh;
	return 0;
} 
