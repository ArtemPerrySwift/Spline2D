#include "Spline2D.h"

std::istream& InterpFuncData::readData(std::istream& in)
{
	in >> point >> funcValue >> weightCoeff;

	if (in.fail())
	{
		programlog::writeErr("Unable to read " + getTypeDataName());
	}

	if (weightCoeff < 0)
		weightCoeff = abs(1 / funcValue);

	return in;
}

std::string InterpFuncData::getTypeDataName()
{
	return "interpolation function data";
}

int CubeErmitBasicFunct2D::u(int i)
{
	i++;
	return 2 * ((i - 1) / 4 % 2) + ((i - 1) % 2);
}

int CubeErmitBasicFunct2D::v(int i)
{
	i++;
	return 2 * ((i - 1) / 8) + (((i - 1) / 2 % 2));
}

double CubeErmitBasicFunct1D::cubeErmitBasicFunct1D(double ksi, int iFunct)
{
	switch (iFunct)
	{
	case 0:
		return 1 - 3 * ksi * ksi + 2 * ksi * ksi * ksi;
	case 1:
		return ksi - 3 * ksi * ksi + ksi * ksi * ksi;
	case 2:
		return 3 * ksi * ksi - 2 * ksi * ksi * ksi;
	case 3:
		return -ksi * ksi + ksi * ksi * ksi;
	default:
		programlog::writeErr("The number of cube basic function should be in [0, 4]");
		return 0;
	}
}

void CubeErmitBasicFunct1D::getStiffMatrix(double hKsi, double G[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D])
{
	G[0][0] = 36/(30* hKsi);

	G[0][1] = G[1][0] = 3;
	G[1][1] = 4 * hKsi;

	G[0][2] = G[2][0] = -36 / hKsi;
	G[1][2] = G[2][1] = -3.0/30;
	G[2][2] = 36/(30.0*hKsi);

	G[0][3] = G[3][0] = 3.0 / 30;
	G[1][3] = G[3][1] = -hKsi / 30;
	G[2][3] = G[3][2] = -3.0 / 30;
	G[3][3] = 4 * hKsi / 30;
}

void CubeErmitBasicFunct1D::getMassMatrix(double hKsi, double M[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D])
{
	M[0][0] = 156 * hKsi / 420;

	M[0][1] = M[1][0] = 22 * hKsi* hKsi/420;
	M[1][1] = 4 * hKsi* hKsi * hKsi / 420;

	M[0][2] = M[2][0] = 54 * hKsi / 420;
	M[1][2] = M[2][1] = 13 * hKsi * hKsi / 420;
	M[2][2] = 156 * hKsi / 420;

	M[0][3] = M[3][0] = -13 * hKsi * hKsi / 420;
	M[1][3] = M[3][1] = -3 * hKsi * hKsi * hKsi / 420;
	M[2][3] = M[3][2] = -22 * hKsi * hKsi / 420;
	M[3][3] = 4 * hKsi * hKsi * hKsi / 420;
}

void CubeErmitBasicFunct1D::getSecondMatrix(double hKsi, double S[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D])
{
	S[0][0] = 60/(hKsi * hKsi * hKsi);

	S[0][1] = S[1][0] = 30/(hKsi * hKsi);
	S[1][1] = 16/hKsi;

	S[0][2] = S[2][0] = -60 / (hKsi * hKsi * hKsi);
	S[1][2] = S[2][1] = -30 / (hKsi * hKsi);
	S[2][2] = 60 / (hKsi * hKsi * hKsi);

	S[0][3] = S[3][0] = 30 / (hKsi * hKsi);
	S[1][3] = S[3][1] = 14 / hKsi;
	S[2][3] = S[3][2] = -30 / (hKsi * hKsi);
	S[3][3] = 16 / hKsi;
}

void CubeErmitBasicFunct2D::getStiffMatrix2D(double hKsi, double hNu, double G[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D])
{
	double Gksi[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];
	double Gnu[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];
	double Mksi[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];
	double Mnu[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];

	getStiffMatrix(hKsi, Gksi);
	getStiffMatrix(hNu, Gnu);
	getMassMatrix(hKsi, Mksi);
	getMassMatrix(hNu, Mnu);

	for (int i = 0; i < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; i++)
	{
		for (int j = i; j < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; j++)
		{
			G[i][j] = G[j][i] = Gksi[u(i)][u(j)] * Mnu[v(i)][v(j)] + Mksi[u(i)][u(j)] * Gnu[v(i)][v(j)];
		}
	}
}

void CubeErmitBasicFunct2D::getMassMatrix2D(double hKsi, double hNu, double M[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D])
{
	double Mksi[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];
	double Mnu[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];

	getMassMatrix(hKsi, Mksi);
	getMassMatrix(hNu, Mnu);

	for (int i = 0; i < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; i++)
	{
		for (int j = i; j < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; j++)
		{
			M[i][j] = M[j][i] = Mksi[u(i)][u(j)] * Mnu[v(i)][v(j)];
		}
	}
}

void CubeErmitBasicFunct2D::getSecondMatrix2D(double hKsi, double hNu, double S[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D])
{
	double Sksi[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];
	double Snu[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];
	double Mksi[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];
	double Mnu[NUM_CUBE_ERMIT_BASIC_FUNCT_1D][NUM_CUBE_ERMIT_BASIC_FUNCT_1D];

	getSecondMatrix(hKsi, Sksi);
	getSecondMatrix(hNu, Snu);
	getMassMatrix(hKsi, Mksi);
	getMassMatrix(hNu, Mnu);

	for (int i = 0; i < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; i++)
	{
		for (int j = i; j < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; j++)
		{
			S[i][j] = S[j][i] = Sksi[u(i)][u(j)] * Mnu[v(i)][v(j)] + Mksi[u(i)][u(j)] * Snu[v(i)][v(j)];
		}
	}
}

double CubeErmitBasicFunct2D::cubeErmitBasicFunct2D(double ksi, double nu, int iFunct)
{
	return cubeErmitBasicFunct1D(ksi, u(iFunct)) * cubeErmitBasicFunct1D(nu, v(iFunct));
}