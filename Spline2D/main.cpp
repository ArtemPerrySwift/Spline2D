#include "Spline2D.h"
#include "DataGenerator.h"

int main()
{
	DataGenerator dataGenerator;
	RegularMesh regularDataGeneratorMesh("InterpFunctMesh.txt");
	dataGenerator.writePureDataInFile(regularDataGeneratorMesh.vertices, "InterpFunctData.txt");

	return 0;
}