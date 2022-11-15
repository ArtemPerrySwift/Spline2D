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

double CubeErmitBasicFunct1D::ksi(double x, double xi, double hx)
{
	return (x - xi) / hx;
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

void CubeErmitBasicFunct2D::getSecondMatrix2D(double hKsi, double hNu, double S[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D])
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

void Spline2D::distrInterpFuncDataToFinElems()
{
	int numInterpFuncData = interpFuncData_s.size();
	interFuncDataOfFinElem.resize(regularFinitMesh.finitElements.size());
	for (int i = 0; i < numInterpFuncData; i++)
	{
		int iFinElemWithInterpFuncData = regularFinitMesh.getFinElemNumByPoint(interpFuncData_s[i].point);
		if (iFinElemWithInterpFuncData > 0)
			interFuncDataOfFinElem[iFinElemWithInterpFuncData].push_back(i);
	}
}

void Spline2D::getLocalA(double A[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem)
{
	double xl, xr;
	double yl, yr;
	FinElem finitElement = regularFinitMesh.finitElements[iFinElem];
	xl = regularFinitMesh.vertices[finitElement.verInd[0]].x;
	xr = regularFinitMesh.vertices[finitElement.verInd[1]].x;
	yl = regularFinitMesh.vertices[finitElement.verInd[0]].y;
	yr = regularFinitMesh.vertices[finitElement.verInd[3]].y;
	double hx = xr - xl;
	double hy = yr - yl;
	int iInterpFincData;
	Coord2D interFuncDataPoint;
	double ksi, nu;
	for (int i = 0; i < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; i++)
	{
		for (int j = 0; j < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; j++)
		{
			for (int k = 0; k < interFuncDataOfFinElem[iFinElem].size(); k++)
			{
				iInterpFincData = interFuncDataOfFinElem[iFinElem][k];
				interFuncDataPoint = interpFuncData_s[iInterpFincData].point;
				ksi = cubeErmitBasicFunct2D.ksi(interFuncDataPoint.x, xl, hx);
				nu = cubeErmitBasicFunct2D.ksi(interFuncDataPoint.y, yl, hy);
				A[i][j] = interpFuncData_s[iInterpFincData].weightCoeff * cubeErmitBasicFunct2D.cubeErmitBasicFunct2D(ksi, nu, i) * cubeErmitBasicFunct2D.cubeErmitBasicFunct2D(ksi, nu, j);
			}
		}
	}
	
}

void Spline2D::getLocalB(double B[NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem)
{
	double xl, xr;
	double yl, yr;
	FinElem finitElement = regularFinitMesh.finitElements[iFinElem];
	xl = regularFinitMesh.vertices[finitElement.verInd[0]].x;
	xr = regularFinitMesh.vertices[finitElement.verInd[1]].x;
	yl = regularFinitMesh.vertices[finitElement.verInd[0]].y;
	yr = regularFinitMesh.vertices[finitElement.verInd[3]].y;
	double hx = xr - xl;
	double hy = yr - yl;
	int iInterpFincData;
	Coord2D interFuncDataPoint;
	double ksi, nu;
	for (int i = 0; i < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; i++)
	{
		for (int k = 0; k < interFuncDataOfFinElem[iFinElem].size(); k++)
		{
			iInterpFincData = interFuncDataOfFinElem[iFinElem][k];
			interFuncDataPoint = interpFuncData_s[iInterpFincData].point;
			ksi = cubeErmitBasicFunct2D.ksi(interFuncDataPoint.x, xl, hx);
			nu = cubeErmitBasicFunct2D.ksi(interFuncDataPoint.y, yl, hy);
			B[i] = interpFuncData_s[iInterpFincData].weightCoeff * cubeErmitBasicFunct2D.cubeErmitBasicFunct2D(ksi, nu, i) * interpFuncData_s[iInterpFincData].funcValue;
		}
	}
}

void Spline2D::getLocalG(double G[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem)
{
	double xl, xr;
	double yl, yr;
	FinElem finitElement = regularFinitMesh.finitElements[iFinElem];
	xl = regularFinitMesh.vertices[finitElement.verInd[0]].x;
	xr = regularFinitMesh.vertices[finitElement.verInd[1]].x;
	yl = regularFinitMesh.vertices[finitElement.verInd[0]].y;
	yr = regularFinitMesh.vertices[finitElement.verInd[3]].y;
	double hx = xr - xl;
	double hy = yr - yl;
	cubeErmitBasicFunct2D.getStiffMatrix2D(hx, hy, G);
	
}

void Spline2D::getLocalS(double S[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem)
{
	double xl, xr;
	double yl, yr;
	FinElem finitElement = regularFinitMesh.finitElements[iFinElem];
	xl = regularFinitMesh.vertices[finitElement.verInd[0]].x;
	xr = regularFinitMesh.vertices[finitElement.verInd[1]].x;
	yl = regularFinitMesh.vertices[finitElement.verInd[0]].y;
	yr = regularFinitMesh.vertices[finitElement.verInd[3]].y;
	double hx = xr - xl;
	double hy = yr - yl;
	cubeErmitBasicFunct2D.getSecondMatrix2D(hx, hy, S);
}

void Spline2D::getGlobalInd(int globalInd[NUM_CUBE_ERMIT_BASIC_FUNCT_2D], int iFinElem)
{
	FinElem finitElement = regularFinitMesh.finitElements[iFinElem];
	for (int i = 0, j = 0; i < QUAD_VER; i++)
	{
		int curInd = 4*finitElement.verInd[i];
		for (int k = 0; k < 4; k++, j++)
			globalInd[j] = curInd + k;
	}
}

void Spline2D::assembleGlobalMatrix()
{
	int ktr = regularFinitMesh.finitElements.size();
	int materialInd;
	//Material material;
	slae.A.fillMatrix(0);

	char x = -37;
	int del = ktr / 35;
	double localG[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D];
	double localA[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D];
	double localS[NUM_CUBE_ERMIT_BASIC_FUNCT_2D][NUM_CUBE_ERMIT_BASIC_FUNCT_2D];
	double localB[NUM_CUBE_ERMIT_BASIC_FUNCT_2D];
	int globalInd[NUM_CUBE_ERMIT_BASIC_FUNCT_2D];
	//std::cout << del;
	//std::cout << setlocale(LC_ALL, NULL) << std::endl;
	//setlocale(LC_ALL, "C");
	//std::cout << setlocale(LC_ALL, NULL) << std::endl;
	//std::cout << std::endl;
	//programlog::ProgressBar progressBar;
	//progressBar.begin();
	for (int i = 0; i < ktr; i++)
	{
		//progressBar.showProgress(double(i + 1) / ktr * 100.0);
		int j, k;
		//materialInd = mesh.finalElementMaterialStorage.finalElemMaterials[i];
		//material = mesh.materialStorage.findMaterialByNum(materialInd);
		//std::cout << "\t Построение локальной матрицы жёсткости" << std::endl;
		//finalElemBase.buildLocalG(localG);
		getGlobalInd(globalInd, i);
		getLocalA(localA, i);
		getLocalG(localG, i);
		getLocalS(localS, i);
		getLocalB(localB, i);
		addLocalMatrixesToGlobalOne(localA, localG, localS, localB, globalInd);
		/*
		std::cout << "Local G" << std::endl;
		double sum;
		for (j = 0; j < 8; j++)
		{
			sum = 0;
			for (k = 0; k < 8; k++)
			{
				sum += localG[j][k];
				std::cout << localG[j][k] << " ";
			}
			std::cout << sum << std::endl;
		}
		*/
		//std::cout << "\t Построение локального вектора b" << std::endl;
		//finalElemBase.buildLocalB(localB);
		//finalElemBase.buildLocalM(localM, 0);

		/*
		std::cout << "Local M" << std::endl;
		for (j = 0; j < 8; j++)
		{
			sum = 0;
			for (k = 0; k < 8; k++)
			{
				sum += localM[j][k];
				std::cout << localM[j][k] << " ";
			}
			std::cout << sum << std::endl;
		}
		std::cout << "Local B" << std::endl;
		for (j = 0; j < 8; j++)
		{
			std::cout << localB[j] << " ";
		}
		*/
		//addLocalToGlobalAll(slae, localG, localM, localB, finalElemBase.globalInd);
	}
	//progressBar.end();
	//std::cout << std::endl;
	//getchar();
}

void Spline2D::addLocalMatrixesToGlobalOne(double A_local[FREE_DEG_2D][FREE_DEG_2D], double G_local[FREE_DEG_2D][FREE_DEG_2D], double S_local[FREE_DEG_2D][FREE_DEG_2D], double B_local[FREE_DEG_2D], int L[FREE_DEG_2D])
{
	int m, l;
	int ind;
	double* di = slae.A.di;
	double* gg = slae.A.gg;
	double* b = slae.b;
	for (m = 0; m < FREE_DEG_2D; m++)
	{
		di[L[m]] += alpha*G_local[m][m] + A_local[m][m] + betta * S_local[m][m];
		b[L[m]] += B_local[m];
	}
	int s, d;
	for (m = 1; m < FREE_DEG_2D; m++)
	{
		for (l = 0; l < m; l++)
		{
			ind = slae.A.getElemIndG(L[m], L[l]);
			gg[ind] += alpha * G_local[m][l] + A_local[m][l] + betta * S_local[m][l];
		}
	}
}

void Spline2D::countSpline()
{
	assembleGlobalMatrix();
	q.resize(slae.A.n);
	los.solve(slae.A, slae.b, q.data(), 1000, 1e-12);
}

double Spline2D::getSplineValue(Coord2D point)
{
	double xl, xr;
	double yl, yr;
	int iFinElem = regularFinitMesh.getFinElemNumByPoint(point);
	if (iFinElem < 0)
		return 0;

	FinElem finitElement = regularFinitMesh.finitElements[iFinElem];
	int globalInd[NUM_CUBE_ERMIT_BASIC_FUNCT_2D];
	getGlobalInd(globalInd, iFinElem);
	xl = regularFinitMesh.vertices[finitElement.verInd[0]].x;
	xr = regularFinitMesh.vertices[finitElement.verInd[1]].x;
	yl = regularFinitMesh.vertices[finitElement.verInd[0]].y;
	yr = regularFinitMesh.vertices[finitElement.verInd[3]].y;
	double hx = xr - xl;
	double hy = yr - yl;
	int iInterpFincData;
	Coord2D interFuncDataPoint;
	double ksi, nu;
	double res = 0;
	for (int i = 0; i < NUM_CUBE_ERMIT_BASIC_FUNCT_2D; i++)
	{
		ksi = cubeErmitBasicFunct2D.ksi(interFuncDataPoint.x, xl, hx);
		nu = cubeErmitBasicFunct2D.ksi(interFuncDataPoint.y, yl, hy);
		res += q[globalInd[i]] * cubeErmitBasicFunct2D.cubeErmitBasicFunct2D(ksi, nu, i);
	}

	return res;
}

Spline2D::Spline2D(std::string fileNameMesh, std::string fileNameFunc, std::string fileNameSpline) : regularFinitMesh(fileNameMesh)
{
	std::ifstream inFunc;
	inFunc.open(fileNameFunc);
	inFunc >> interpFuncData_s;
	inFunc.close();

	std::ifstream inSpline;
	inSpline.open(fileNameSpline);
	inSpline >> alpha >> betta;
	if (inSpline.fail())
		programlog::writeErr("Unable to read alpha and betta coefficients for spline");
	
	if(alpha < 0)
		programlog::writeErr("Alpha coefficient for spline connot be < 0");

	if (betta < 0)
		programlog::writeErr("Betta coefficient for spline connot be < 0");

	inSpline.close();

	countSpline();
}
