#pragma once
#include <array>
#include "InData.h"
#include "Coord2D.h"

const int QUAD_VER = 4;

struct FinElem
{
	std::array<int, QUAD_VER> verInd;
};

struct SepParam : InData
{
	int n;
	double q;

	SepParam();
	SepParam(int n, double q);

	virtual std::istream& readData(std::istream& in) override;
	virtual std::string getTypeDataName() override;
};

struct AxisCoordinates : InData
{
public:
	AxisCoordinates(std::istream& in, std::string axisName);
	std::istream& init(std::istream& in);
	virtual std::string getTypeDataName() override;
	std::string getAxisName();
	double operator[] (int index);
	int count();
private:
	std::vector<double> coordinates;
	std::vector<double> basicCoordinates;
	std::vector<SepParam> sepParams;
	std::string axisName;

	int getNumCoordsOfSeperatedAxis();
	int getCoordsOfSeperatedAxis();
	void clearHelpVectors();
	virtual std::istream& readData(std::istream& in) override;
};

std::istream& operator >> (std::istream& in, AxisCoordinates& axisCoordinates);
/*
struct MeshInData : InData
{
public:
	AxisCoordinates coordAxisXParams;
	AxisCoordinates coordAxisYParams;

	MeshInData(std::string fileName);
	MeshInData(std::istream& in);
	virtual std::istream& readData(std::istream& in) override;
	virtual std::string getTypeDataName() override;

};
*/
struct RegularFinitMesh
{
public:
	std::vector<Coord2D> vertices;
	std::vector<FinElem> finitElements;

	RegularFinitMesh(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
	RegularFinitMesh(std::istream& in);
	RegularFinitMesh(std::string fileName);
	int getFinElemNumByPoint(Coord2D point);
private:
	int nXCoords, nYCoords;
	int nXFinElems, nYFinElems;

	void init(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
	void init(std::istream& in);
	int getGlobalVertNum(int ix, int iy);
	void fillVertices(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
	void fillFinitElements(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
};