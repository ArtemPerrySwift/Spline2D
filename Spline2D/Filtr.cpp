#include "Filtr.h"


void Filtr::filtrate(Spline2D& spline, double errMaxCoeff)
{
	interpFuncData_s = spline.getInterpFunctData();
	do
		maxErr = errMaxCoeff*calcErrMean(spline);
	while (reduceFuncWeights(spline) > 0);
	interpFuncData_s.clear();
}

double Filtr::calcErrMean(Spline2D& spline)
{
	int nData = interpFuncData_s.size();
	double errSum = 0;
	for (int i = 0; i < nData; i++)
	{
		errSum += abs(spline.getSplineValue(interpFuncData_s[i].point) - interpFuncData_s[i].funcValue);
	}
	return errSum / nData;
}

int Filtr::reduceFuncWeights(Spline2D& spline)
{
	int nData = interpFuncData_s.size();
	double err;
	WeightToChange weightToChange;
	for (int i = 0; i < nData; i++)
	{
		err = abs(spline.getSplineValue(interpFuncData_s[i].point) - interpFuncData_s[i].funcValue);
		if (err > maxErr)
		{
			interpFuncData_s[i].weightCoeff /= 2 * err / maxErr;
			weightToChange.indWeight = i;
			weightToChange.newWeightValue = interpFuncData_s[i].weightCoeff;
			weightsToChange.push_back(weightToChange);
		}
	}
	
	spline.changeFuncWeightValues(weightsToChange);
	int changedWeightsNum = weightsToChange.size();
	weightsToChange.clear();
	return changedWeightsNum;
}