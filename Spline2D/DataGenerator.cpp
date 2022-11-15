#include "DataGenerator.h"

double DataGenerator::getPureFunctValue(Coord2D point)
{
	return point.x;
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
	out << pointsNum;
	for (int i = 0; i < pointsNum; i++)
	{
		out << points[i].x << " " << points[i].y << " " << getPureFunctValue(points[i]) << -1;
	}
	out.close();
}

void DataGenerator::writeDirtDataInFile(std::vector<Coord2D>& points, std::string fileName)
{
	std::ofstream out;
	out.open(fileName);
	int pointsNum = points.size();
	out << pointsNum;
	for (int i = 0; i < pointsNum; i++)
	{
		out << points[i].x << " " << points[i].y << " " << getDirtFunctValue(points[i]) << -1;
	}
	out.close();
}
