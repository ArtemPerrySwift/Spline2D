#pragma once
#include <array>
#include "InData.h"
#include "Coord2D.h"

const int QUAD_VER = 4;
const int BASIC_FUNC_NUM = QUAD_VER * QUAD_VER;

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

struct RegularMesh
{
public:
	std::vector<Coord2D> vertices;
	RegularMesh(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
	RegularMesh(std::istream& in);
	RegularMesh(std::string fileName);
	void writeMeshInFile(std::string fileName);
	void writeMeshInFile(std::ostream& out);
protected:
	int nXCoords, nYCoords;
	int nXFinElems, nYFinElems;

	virtual void init(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
	void init(std::istream& in);
	void fillVertices(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
};

struct RegularFinitMesh : RegularMesh
{
public:
	std::vector<FinElem> finitElements;

	RegularFinitMesh(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
	RegularFinitMesh(std::istream& in);
	RegularFinitMesh(std::string fileName);
	int getFinElemNumByPoint(Coord2D point);
private:
	void init(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates) override;
	void init(std::istream& in);
	int getGlobalVertNum(int ix, int iy);
	void fillFinitElements(AxisCoordinates& AxisXCoordinates, AxisCoordinates& AxisYCoordinates);
};