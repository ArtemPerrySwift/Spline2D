#include "DataGenerator.h"
#include <cmath>

double DataGenerator::getPureFunctValue(Coord2D point)
{
	return std::sin(point.x + point.y);
	//return point.x + point.y;
}

double DataGenerator::getDirtFunctValue(Coord2D point, double dirtLvl)
{
	return getPureFunctValue(point) * (1 + rand()/RAND_MAX*dirtLvl);
}

void DataGenerator::writePureDataInFile(std::vector<Coord2D>& points, std::string fileName)
{
	std::ofstream out;
	out.open(fileName);
	int pointsNum = points.size();
	out << pointsNum << std::endl;
	for (int i = 0; i < pointsNum; i++)
	{
		out << points[i].x << "\t" << points[i].y << "\t" << getPureFunctValue(points[i]) << "\t" << -1 << std::endl;
	}
	out.close();
}

void DataGenerator::writeDirtDataInFile(std::vector<Coord2D>& points, std::string fileName, double dirtLvl)
{
	std::ofstream out;
	out.open(fileName);
	int pointsNum = points.size();
	out << pointsNum << std::endl;
	for (int i = 0; i < pointsNum; i++)
	{
		out << points[i].x << "\t" << points[i].y << "\t" << getDirtFunctValue(points[i], dirtLvl) << "\t" << -1 << std::endl;
	}
	out.close();
}
