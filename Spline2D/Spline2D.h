#pragma once
#include "RegularFinitMesh.h"
#include "slae.h"
#include "los.h"

const int NUM_CUBE_ERMIT_BASIC_FUNCT_1D = 4;
const int NUM_CUBE_ERMIT_BASIC_FUNCT_2D = NUM_CUBE_ERMIT_BASIC_FUNCT_1D * NUM_CUBE_ERMIT_BASIC_FUNCT_1D;

const int FREE_DEG_1D = NUM_CUBE_ERMIT_BASIC_FUNCT_1D;
const int FREE_DEG_2D = NUM_CUBE_ERMIT_BASIC_FUNCT_2D;

struct InterpFuncData : InData
{
	Coord2D point;
	double funcValue;
	double weightCoeff;

	virtual std::istream& readData(std::istream& in) override;
	virtual std::string getTypeDataName() override;
};

struct WeightToChange
{
	unsigned int indWeight;
	double newWeightValue;
};

class CubeErmitBasicFunct1D
{
public:
	double ksi(double x, double xi, double hx);
	double cubeErmitBasicFunct1D(double ksi, int iFunct);
	double cubeErmitBasicDifFunct1D(double ksi, int iFunct);
	void getStiffMatrix(double hKsi, double G[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D]);
	void getMassMatrix(double hKsi, double M[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D]);
	void getSecondMatrix(double hKsi, double S[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D]);
};

class CubeErmitBasicFunct2D : public CubeErmitBasicFunct1D
{
public:
	double cubeErmitBasicFunct2D(double ksi, double nu, int iFunct);
	double cubeErmitBasicDifKsiFunct2D(double ksi, double nu, int iFunct);
	double cubeErmitBasicDifNuFunct2D(double ksi, double nu, int iFunct);
	void getStiffMatrix2D(double hKsi, double hNu, double G[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D]);
	void getMassMatrix2D(double hKsi, double hNu, double M[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D]);
	void getSecondMatrix2D(double hKsi, double hNu, double S[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D]);

private:
	int u(int i);
	int v(int i);
};

class Spline2D
{
public:
	friend class SplineCalculator;

	Spline2D(std::string fileNameMesh, std::string fileNameFunc, std::string fileNameSpline);
	Spline2D(std::string fileNameMesh, std::string fileNameSpline, std::vector<InterpFuncData> interpFuncData_s);

	double getSplineValue(Coord2D point);
	double getDifXSplineValue(Coord2D point);
	double getDifYSplineValue(Coord2D point);
	void writeSplineValuesInFile(std::vector<Coord2D> points, std::string fileName);
	void writeDifSplineValuesInFile(std::vector<Coord2D> points, std::string fileName);
	std::vector<InterpFuncData> getInterpFunctData();
	void changeFuncWeightValue(WeightToChange weightToChange);
	void changeFuncWeightValues(std::vector<WeightToChange>& weightsToChange);

private:
	RegularFinitMesh regularFinitMesh;
	CubeErmitBasicFunct2D cubeErmitBasicFunct2D;
	slae::SLAE<SparseMatrixSym> slae;
	std::vector<InterpFuncData> interpFuncData_s;
	std::vector<std::vector<int>> interFuncDataOfFinElem;
	std::vector<double> q;
	LOS los;
	double alpha;
	double betta;

	void getLocalA(double A[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem);
	void getLocalB(double B[NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem);
	void getLocalG(double G[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem);
	void getLocalS(double S[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem);
	void getGlobalInd(int globalInd[NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem);
	void distrInterpFuncDataToFinElems();
	void addLocalMatrixesToGlobalOne(double A_local[FREE_DEG_2D][FREE_DEG_2D], double G_local[FREE_DEG_2D][FREE_DEG_2D], double S_local[FREE_DEG_2D][FREE_DEG_2D], double B_local[FREE_DEG_2D], int L[FREE_DEG_2D]);
	void assembleGlobalMatrix();
	void countSpline();
	void readInterpFuncData(std::string fileNameFunc);
	void readSplineParams(std::string fileNameSpline);
	void slaeInit();
	void privChangeFuncWeightValue(WeightToChange weightToChange);
};
