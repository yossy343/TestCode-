// Test_Kriging.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// Wisarut Comments //
#include "pch.h"

#include <iomanip>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <time.h>
#include <ctime>
#include <sys/timeb.h>
#include <chrono>

#include <queue>
#include <functional>

#include <random>
using namespace std;

struct Point
{
	int x;
	int y;
	int z;
};
double Cald(double sill, double nugget, double range, double A)
{
	double result;
	result = nugget + sill * (1 - exp(-3 * A / range));
	return result;
}
void print(vector< vector<double> > A) {
	int n = A.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 2 * n; j++) {
			cout << A[i][j] << "\t";
			if (j == n - 1) {
				cout << "| ";
			}
		}
		cout << "\n";
	}
	cout << endl;
}

void calculateInverse(vector< vector<double>> &A, int size)
{
	int n = size;
	printf("A.size() = %d\n", A.size());
	for (int i = 0; i < n; i++)
	{
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k < n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = A[k][i];
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k < 2 * n; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j < 2 * n; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	for (int i = n - 1; i >= 0; i--) {
		for (int k = n; k < 2 * n; k++) {
			A[i][k] /= A[i][i];
		}
		// this is not necessary, but the output looks nicer:
		A[i][i] = 1;

		for (int rowModify = i - 1; rowModify >= 0; rowModify--) {
			for (int columModify = n; columModify < 2 * n; columModify++) {
				A[rowModify][columModify] -= A[i][columModify]
					* A[rowModify][i];
			}
			// this is not necessary, but the output looks nicer:
			A[rowModify][i] = 0;
		}
	}
}
bool Multiply_Matrix(vector<vector<double>> &m, vector<vector<double>> &n, vector<vector<double>> &Out, int r1, int c1, int r2, int c2)
{
	// Multiplying matrix a and b and storing in array mult.
	for (int i = 0; i < r1; ++i)
	{
		for (int j = 0; j < c2; ++j)
		{
			for (int k = 0; k < c1; ++k)
			{
				Out[i][j] += m[i][k] * n[k][j];
			}
		}
	}
	return true;
}

int main()
{
	int nugget = 1;
	int sill = 4;
	int range = 6;

	printf("Test\n");
	Point P1, P2, P3, P4, P5, P6;
	P1.x = 0; P1.y = 0; P1.z = 1;
	P2.x = 1; P2.y = 0; P2.z = 2;
	P3.x = 2; P3.y = 0; P3.z = 4;
	P4.x = 0; P4.y = 1; P4.z = 5;
	P5.x = 0; P5.y = 2; P5.z = 6;
	P6.x = 2; P6.y = 2; P6.z = 27;

	vector<Point> Neigbore;
	Neigbore.clear();
	Neigbore.push_back(P1);
	Neigbore.push_back(P2);
	Neigbore.push_back(P3);
	Neigbore.push_back(P4);
	Neigbore.push_back(P5);
	Neigbore.push_back(P6);

	int Size = Neigbore.size() + 1;

	vector<double> Distance0toN;
	Distance0toN.clear();
	Point S0;
	S0.x = 1;
	S0.y = 1;
	S0.z = 0;

	for (int i = 0; i < Neigbore.size(); i++)
	{
		double D = sqrt(pow(S0.x - Neigbore.at(i).x, 2) + pow(S0.y - Neigbore.at(i).y, 2));
		Distance0toN.push_back(D);
		printf("D 0 to Neigbore %d = %.3f\n", i + 1, Distance0toN.at(i));
	}
	//Distance0toN.push_back(1); // Add to make 2 matrix equal
	printf("----------------------------------------------\n");
	// For n neigbore, create matrix with size (n+1)x(n+1)
	printf("Size = %d\n", Size);
	vector<vector<double>> DistanceOfEachPoint(Size, vector<double>(Size, 0));
	vector<vector<double>> Matrix_A(Size, vector<double>(Size * 2, 0));
	vector<vector<double>> Matrix_A0(Size, vector<double>(Size * 2, 0));
	vector<vector<double>> Matrix_AInv(Size, vector<double>(Size, 0));
	vector<vector<double>> Check(Size, vector<double>(Size, 0));

	vector<vector<double>> C(Size, vector<double>(1, 0));
	vector<vector<double>> W(Size, vector<double>(1, 0));

	//DistanceOfEachPoint.clear();
	//Matrix_A.clear();
	//Matrix_A0.clear();
	//Matrix_AInv.clear();
	//Check.clear();
	//C.clear();
	//W.clear();

	//printf("Size = %d, %d \n", Matrix_A.size(), Matrix_A.at(0).size());
	int i, j;
	printf("Matrix Distance\n");
	for (i = 0; i < Size; i++)
	{
		for (j = 0; j < Size; j++)
		{
			printf("[%d][%d]\t",i,j);
			if (i < Size - 1 && j < Size - 1)
				DistanceOfEachPoint[i][j] = sqrt(pow(Neigbore[i].x - Neigbore[j].x, 2) + pow(Neigbore[i].y - Neigbore[j].y, 2));
			else if (i == Size - 1)
				DistanceOfEachPoint[i][j] = 1;
			else if (j == Size - 1)
				DistanceOfEachPoint[i][j] = 1;
		    else
			{

			}
			if (i == Size - 1 && j == Size - 1)
			{
				DistanceOfEachPoint[i][j] = 0;
			}
			printf("%.3f\t", DistanceOfEachPoint[i][j]);
			
		}
		printf("\n");
	}

//	Sleep(1000000);
	printf("Matrix Distance Done\n");
	printf("----------------------------------------------\n");


	printf("Check point\n");
	i = 0;
	j = 0;
	printf("Matrix Gamma\n");

	for (i = 0; i < Size; i++)
	{
		for (j = 0; j < Size; j++)
		{
			//printf("[%d][%d]\t",i,j);
			if (i < Size - 1 && j < Size - 1)
			{
				Matrix_A[i][j] = Cald(sill, nugget, range, DistanceOfEachPoint[i][j]);
				Matrix_A0[i][j] = Matrix_A[i][j];
			}

			else if (i == Size - 1)
			{
				Matrix_A[i][j] = 1;
				Matrix_A0[i][j] = Matrix_A[i][j];
			}
			else if (j == Size - 1)
			{
				Matrix_A[i][j] = 1;
				Matrix_A0[i][j] = Matrix_A[i][j];
			}

			else
			{

			}
			//
			if (i == j)
			{
				Matrix_A[i][j] = 0;
				Matrix_A0[i][j] = 0;
			}
			printf("%.3f\t", Matrix_A[i][j]);
			//printf("[%d][%d]\t",i,j);
		}
		printf("\n");
	}
	printf("----------------------------------------------\n");
	int n = Size;
	for (int i = 0; i < n; i++)
	{
		Matrix_A[i][n + i] = 1;

	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2 * n; j++)
		{
			printf("%.2f\t", Matrix_A[i][j]);
		}
		printf("\n");
	}
	printf("----------------------------------------------\n");
	calculateInverse(Matrix_A, Size);
	printf("Matrix A Inverse\n");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Matrix_AInv[i][j] = Matrix_A[i][j + n];
			printf("%.2f\t", Matrix_AInv[i][j]);
		}
		printf("\n");
	}
	printf("----------------------------------------------\n");
	printf("Check Inv Value\n");
	Multiply_Matrix(Matrix_A0, Matrix_AInv, Check, 7, 7, 7, 7);
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%.2f\t", Check[i][j]);
		}
		printf("\n");
	}
	printf("----------------------------------------------\n");
	
	printf("D P0 to P[i]\n");
	printf("Size of Matrix C = %d\n",C.size());
	printf("Size of Matrix Distance0toN = %d\n", Distance0toN.size());
	for (i = 0; i < Size; i++)
	{
		if (i < Size - 1)
		{
			C[i][0] = Cald(sill, nugget, range, Distance0toN[i]);
			printf("Index = %d\tDistance0toN = %.3f\tC[0 to N][%d] = %.3f \n", i, Distance0toN[i], i, C[i][0]);
		}    	
	}
	C[Size - 1][0] = 1;

	printf("----------------------------------------------\n");
	printf("Weight of P[i]\n");

	Multiply_Matrix(Matrix_AInv, C, W, 7, 7, 7, 1);

	printf("Size of Matrix C = %d\n", C.size());
	printf("Size of Matrix W = %d\n\n", W.size());
	for (i = 0; i < Size; i++)
	{
		//W[i][0] = W[i][0] - 1;
		printf("W[0 to N][%d] = %.3f \n", i, W[i][0]);
	}
	printf("----------------------------------------------\n");
	printf("Predicted Value\n");
	double PredictVal = 0;

	for (i = 0; i < Size-1; i++)
	{
		PredictVal += W[i][0] * Neigbore[i].z;

	}
	printf("PredictVal = %.3f \n", PredictVal);
	return 0;
}
