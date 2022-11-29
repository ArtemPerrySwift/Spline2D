#include "Spline2D.h"
#include "DataGenerator.h"

int main()
{
	DataGenerator dataGenerator;
	RegularMesh regularDataGeneratorMesh("InterpFunctMesh.txt");
	dataGenerator.writePureDataInFile(regularDataGeneratorMesh.vertices, "InterpFunctData.txt");
	Spline2D spline("SplineMesh.txt", "InterpFunctData.txt", "SplineParameters.txt");

	RegularMesh regularSolutionToFindMesh("ToFindSolutionMesh.txt");
	spline.writeSplineValuesInFile(regularSolutionToFindMesh.vertices, "SplineValues.txt");
	spline.writeDifSplineValuesInFile(regularSolutionToFindMesh.vertices, "SplineDifValues.txt");

	return 0;
}