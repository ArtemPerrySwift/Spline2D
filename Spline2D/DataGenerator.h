#pragma once
#include "RegularFinitMesh.h"
#include <random>
class DataGenerator
{
public:
	double getPureFunctValue(Coord2D point);
	double getDirtFunctValue(Coord2D point, double dirtLvl);
	void writePureDataInFile(std::vector<Coord2D> points, std::string fileName);
	void writeDirtDataInFile(std::vector<Coord2D> points, std::string fileName);
	
};
