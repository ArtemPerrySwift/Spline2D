#pragma once
#include "RegularFinitMesh.h"
#include <functional>

const int NUM_CUBE_ERMIT_BASIC_FUNCT_1D = 4;
const int NUM_CUBE_ERMIT_BASIC_FUNCT_2D = NUM_CUBE_ERMIT_BASIC_FUNCT_1D * NUM_CUBE_ERMIT_BASIC_FUNCT_1D;

struct InterpFuncData : InData
{
	Coord2D point;
	double funcValue;
	double weightCoeff;

	virtual std::istream& readData(std::istream& in) override;
	virtual std::string getTypeDataName() override;
};

class CubeErmitBasicFunct1D
{
public:
	double cubeErmitBasicFunct1D(double ksi, int iFunct);
	void getStiffMatrix(double hKsi, double G[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D]);
	void getMassMatrix(double hKsi, double M[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D]);
	void getSecondMatrix(double hKsi, double S[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D]);
};

class CubeErmitBasicFunct2D : public CubeErmitBasicFunct1D
{
public:
	double cubeErmitBasicFunct2D(double ksi, double nu, int iFunct);
	void getStiffMatrix2D(double hKsi, double hNu, double G[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D]);
	void getMassMatrix2D(double hKsi, double hNu, double M[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D]);
	void getSecondMatrix2D(double hKsi, double hNu, double S[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D]);

private:
	int u(int i);
	int v(int i);
};

class Spline2D
{
	RegularFinitMesh regularFinitMesh;
	CubeErmitBasicFunct2D cubeErmitBasicFunct2D;
};
