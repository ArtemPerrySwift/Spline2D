#include "RegularFinitMesh.h"

SepParam::SepParam()
{
	n = 0;
	q = 0.0;
}

SepParam::SepParam(int n, double q)
{
	this->n = n;
	this->q = q;
}

std::istream& SepParam::readData(std::istream& in)
{
	in >> n >> q;

	if (in.fail())
	{
		programlog::writeErr("Unable to read " + getTypeDataName());
	}

	if (n < 1)
	{
		programlog::writeErr("Number of intervals in seperation parameters cannot be < 1");
	}

	if (q < 0)
	{
		q = -1 / q;
	}

	return in;
}

std::string SepParam::getTypeDataName()
{
	return "seperation parameters";
}


std::string AxisCoordinates::getTypeDataName()
{
	return "coordinate axis parameters";
}

std::istream& AxisCoordinates::readData(std::istream& in)
{
	int nCoords;
	in >> nCoords;

	if (in.fail())
		programlog::writeErr("Unable to read number of coordinates of " + axisName + " axis");
	

	if (nCoords < 2)
		programlog::writeErr("Number of coordinates of " + axisName + " axis cannot be < 2");
	

	int nIntervals = nCoords - 1;

	basicCoordinates.resize(nCoords);
	sepParams.resize(nIntervals);

	in >> basicCoordinates[0];
	if (in.fail())
		programlog::writeErr("Unable to read begin coordinate of " + axisName + " axis");
	

	for (int i = 0; i < nIntervals; i++)
	{
		in >> sepParams[i];
		in >> basicCoordinates[i + 1];
		if (in.fail())
			programlog::writeErr("Unable to read "+ std::to_string(i + 1) + " coordinate of " + axisName + " axis");
		if(basicCoordinates[i] < basicCoordinates[i + 1])
			programlog::writeErr("Error of " + std::to_string(i + 1) + " coordinate of " + axisName + " axis : Next coordinate cannot be less then previous one");
	}

	return in;
}

std::istream& AxisCoordinates::init(std::istream& in)
{
	std::istream& in = readData(in);
	getCoordsOfSeperatedAxis();
	clearHelpVectors();
	return in;
}

std::istream& operator >> (std::istream& in, AxisCoordinates& axisCoordinates)
{
	return axisCoordinates.init(in);
}

AxisCoordinates::AxisCoordinates(std::istream& in, std::string axisName)
{
	this->axisName = axisName;
	readData(in);
}

std::string AxisCoordinates::getAxisName()
{
	return axisName;
}

int AxisCoordinates::getNumCoordsOfSeperatedAxis()
{
	int sepParamsNum = sepParams.size();
	int numCoordsOfSeperatedAxis = 0;
	for (int i = 0; i < sepParamsNum; numCoordsOfSeperatedAxis += sepParams[i].n);
	return numCoordsOfSeperatedAxis;
}

int AxisCoordinates::getCoordsOfSeperatedAxis()
{
	coordinates.resize(getNumCoordsOfSeperatedAxis());
	int nIntervals = sepParams.size();
	int coordInd;
	for (int i = 0, coordInd = 0; i < nIntervals; i++)
	{
		int nSubIntervals = sepParams[i].n;
		int dischRatio = sepParams[i].q;
		double step;
		if (dischRatio == 1.0)
		{
			step = (basicCoordinates[i + 1] - basicCoordinates[i]) / dischRatio;
			for (int j = 0; j < nSubIntervals; j++, coordInd++)
				coordinates[coordInd] = basicCoordinates[i] + i * step;		
		}
		else
		{
			double bufCoeff = (1 - dischRatio) / (1 - pow(dischRatio, nSubIntervals));
			double step = (basicCoordinates[i + 1] - basicCoordinates[i]) * bufCoeff;

			coordinates[coordInd] = basicCoordinates[i];
			coordInd++;

			for (int j = 1; j < nSubIntervals; j++, coordInd++, step *= dischRatio)
				coordinates[coordInd] = coordinates[coordInd - 1] + step;
		}

		coordinates[coordInd] = basicCoordinates[i + 1];
		coordInd++;
	}
}

double AxisCoordinates::operator[] (int index)
{
	return coordinates[index];
}

int AxisCoordinates::count()
{
	return coordinates.size();
}

void AxisCoordinates::clearHelpVectors()
{
	basicCoordinates.clear();
	sepParams.clear();
}

/*
std::istream& MeshInData::readData(std::istream& in)
{
	in >> coordAxisXParams;
	in >> coordAxisYParams;
	return in;
}



MeshInData::MeshInData(std::istream& in) : coordAxisXParams(in, "X"), coordAxisYParams(in, "Y")
{
	readData(in);
}

std::string MeshInData::getTypeDataName()
{
	return "mesh input data";
}
*/

int RegularFinitMesh::getGlobalVertNum(int ix, int iy)
{
	return nXCoords * iy + ix;
}

RegularFinitMesh::RegularFinitMesh(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates)
{
	fillVertices(AxisXCoordinates, AxisYCoordinates);
	fillFinitElements(AxisXCoordinates, AxisYCoordinates);
}

void RegularFinitMesh::fillVertices(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates)
{
	nXCoords = AxisXCoordinates.count();
	nYCoords = AxisYCoordinates.count();

	int nVertices = nXCoords * nYCoords;
	vertices.resize(nVertices);
	Coord2D vertex;
	for (int iy = 0, iVertex = 0; iy < nYCoords; iy++)
	{
		vertex.y = AxisYCoordinates[iy];
		for (int ix = 0; ix < nXCoords; ix++, iVertex++)
		{
			vertex.x = AxisXCoordinates[ix];
			vertices[iVertex] = vertex;
		}
	}
}

void RegularFinitMesh::fillFinitElements(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates)
{
	nXFinElems = AxisXCoordinates.count() - 1;
	nYFinElems = AxisYCoordinates.count() - 1;

	FinElem finitElement;
	int indVertInCurrentYRow = 0;
	int indVertInNextYRow = nXCoords;

	for (int iy = 0, iFinitElement = 0; iy < nXFinElems; iy++)
	{
		for (int ix = 0; ix < nYFinElems; ix++, iFinitElement++, indVertInCurrentYRow++, indVertInNextYRow++)
		{
			finitElement.verInd[0] = indVertInCurrentYRow;
			finitElement.verInd[1] = indVertInCurrentYRow + 1;

			finitElement.verInd[0] = indVertInNextYRow;
			finitElement.verInd[1] = indVertInNextYRow + 1;

			finitElements[iFinitElement] = finitElement;
		}

		indVertInCurrentYRow++;
		indVertInNextYRow++;
	}
}