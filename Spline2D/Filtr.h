#pragma once
#include "Spline2D.h"

class Filtr
{
public:
	void filtrate(Spline2D& spline, double errMaxCoeff);
private:
	double calcErrMean(Spline2D& spline);
	int reduceFuncWeights(Spline2D& spline);
	std::vector<InterpFuncData> interpFuncData_s;
	std::vector<WeightToChange> weightsToChange;
	double maxErr;
};
