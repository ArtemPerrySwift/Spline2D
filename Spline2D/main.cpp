#include "Spline2D.h"
#include "DataGenerator.h"

int main()
{
	DataGenerator dataGenerator;
	RegularMesh regularDataGeneratorMesh("InterpFunctMesh.txt");
	dataGenerator.writePureDataInFile(regularDataGeneratorMesh.vertices, "InterpFunctData.txt");
	Spline2D spline("SplineMesh.txt", "InterpFunctData.txt", "SplineParameters.txt");
	spline.writeSplineValuesInFile(regularDataGeneratorMesh.vertices, "SplineValues.txt");
	return 0;
}