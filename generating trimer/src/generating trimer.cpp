//============================================================================
// Name        : generating.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cmath>


using namespace std;

void matvectmultiply(double mat1[][3], double vect[], double *result)
{
	for(int q = 0; q < 3; q++)
	{
		result[q] = 0.0;
	}
	for (int x = 0; x < 3; x++)
	{
		for (int c = 0; c < 3; c++)
		{
			result[x] += mat1[x][c] * vect[c];
		}
	}
}

double rand_DoubleRange(double a, double b)
{
	return ((b - a)*((double)rand() / RAND_MAX)) + a;
}

int rand_IntRange(int a, int b)
{
	return rand() % (b - a + 1) + a;
}

void randomonsphere(double c[3], double r)
{
    double pi = 3.14159265f;
    double theta = ((double)rand() / (double)RAND_MAX) * 2 * pi;
    double u = 2 * ((double)rand() / (double)RAND_MAX) - 1;
    c[0] = r*cos(theta)*sqrt(1.0 - u*u);
    c[1] = r*sin(theta)*sqrt(1.0 - u*u);
    c[2] = r*u ;

}


int main()
{
	int seed;
	cout << "enter seed (int): "<< endl;
	cin >> seed;
	int number;
	srand(seed);
	double initial[9]= {0.0f,0.0f,0.0f,3.2094f,2.26818f,0.0f,6.4188f,0.0f,0.0f};
	double temprodvect1[3] = { 0, 0, 0 };
	double temprodvect2[3] = { 0, 0, 0 };
	double temprodvect3[3] = { 0, 0, 0 };
	double coordtemp[3];
	cout<< "generating fixed, bent, trimers in random orientations "<< endl;
	cout<< "number of chains? " << endl;
	cin >> number;
	double Pi = 3.14159;
	double theta;
	double phi;
	double gamma;
	double coordlist[9*10000];
	for(int i = 0; i < 10000;i++)
	{
		theta = rand_DoubleRange(-Pi, Pi);
		phi = rand_DoubleRange(-Pi, Pi);
		gamma = rand_DoubleRange(-Pi, Pi);
		double xrm[3][3];
		double yrm[3][3];
		double zrm[3][3];
		xrm[0][0] = 1;
		xrm[0][1] = 0;
		xrm[0][2] = 0;
		xrm[1][0] = 0;
		xrm[1][1] = cos(theta);
		xrm[1][2] = -sin(theta);
		xrm[2][0] = 0;
		xrm[2][1] = sin(theta);
		xrm[2][2] = cos(theta);
		yrm[0][0] = cos(phi);
		yrm[0][1] = 0;
		yrm[0][2] = sin(phi);
		yrm[1][0] = 0;
		yrm[1][1] = 1;
		yrm[1][2] = 0;
		yrm[2][0] = -sin(phi);
		yrm[2][1] = 0;
		yrm[2][2] = cos(phi);
		zrm[0][0] = cos(gamma);
		zrm[0][1] = -sin(gamma);
		zrm[0][2] = 0;
		zrm[1][0] = sin(gamma);
		zrm[1][1] = cos(gamma);
		zrm[1][2] = 0;
		zrm[2][0] = 0;
		zrm[2][1] = 0;
		zrm[2][2] = 1;
		for(int j =0;j<3;j++)
		{
			for(int z =0;z<3;z++)
			{
				coordtemp[z]= initial[j*3+z];
			}

			matvectmultiply(xrm, coordtemp, temprodvect1);
			matvectmultiply(yrm, temprodvect1, temprodvect2);
			matvectmultiply(zrm, temprodvect2, temprodvect3);
			for(int z = 0; z<3;z++)
			{
				coordlist[i*3*3+j*3+z] = temprodvect3[z];
			}
		}
	}
	//for(int i = 0; i < 3*9; i++)
	//{
	//	cout << "index "<< i << " is " << coordlist[i]<<endl;
	//}

	ofstream coordfile("coordinates.dat");
	if (coordfile.is_open())
	{
		for (int i = 0; i < 9*10000; ++i)
		{
			coordfile << coordlist[i]<<endl;
		}
		coordfile.close();
	}
	// vector list
	//generate new array for vectors
	double sited[3];
	int vectlen = 10000;
	double vectorlist[3*10000];
	for (int p = 0; p < vectlen; p++)
	{
	    randomonsphere(sited, 1);
	        for (int c = 0; c < 3; c++)
	        {
	            vectorlist[p*3 + c] = sited[c];
	        }
	}
	//Main 2 Chain Simulation Loop

	//mcmax is number of 2 chain samples
	int mcmax = 200000;

	//the minimum and maximum indexed distances to place chains apart.  dr = 0.1 so this is 3.0 and 36.0 angstroms
	int min = 17;
	int max = 360;
	// other parameters
	int ngrid = 2048;
	double delr = 0.1f;
	// simulation constants
	// T is the temperature in kelvin
	//kb is boltzmann's constant in kCal/mol
	double T = 65.0f;
	double kb = 0.0019833794749f;
	//molecular parameters
	double sigmaA = 3.93f;
	double sigmaB = 3.93f;
	double bl = 3.93f;
	double solFactor = 1.0f;
	int topmol[3] = {0,1,0};
	// makes an array of the types of interactions between two sites on two
	int interactionlist[9];
	double sigmalist[3] = {sigmaA, .5f*(sigmaA + sigmaB),sigmaB};
	//dummy solvation potential
	double wr[3*2048];
	double gr[3*2048];
	for(int i = 0; i< 3*2048; i++)
	{
		wr[i] = 0.0f;
	}
	for(int i = 0; i< 3*2048; i++)
	{
		gr[i] = 0.0f;
	}
	for(int i = 0; i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			if (topmol[i]==0&&topmol[j]==0)
			{
				interactionlist[i+3*j] = 0;
			}
			else if (topmol[i]==0&&topmol[j]==1)
			{
				interactionlist[i+3*j] = 1;
			}
			else if (topmol[i]==1&&topmol[j]==0)
			{
				interactionlist[i+3*j] = 1;
			}
			else if (topmol[i]==1&&topmol[j]==1)
			{
				interactionlist[i+3*j] = 2;
			}
		}
	}
	/*
	for(int i = 0;i<9;i++)
	{
		cout << "element "<< i << " is "<< interactionlist[i]<< endl;
	}
	*/
	double tempsitecoord1[3];
	double tempsitecoord2[3];
	double tempcoordinates1[3*3];
	double tempcoordinates2[3*3];
	double distlist[9];
	int distindexlist[9];
	//main MC loop

	for(int i = 0; i < mcmax; i++)
	{
		for(int n = min;n < max+1; n++)
		{
			int molchoice1 = rand_IntRange(0, 10000-1)*9;
			int molchoice2 = rand_IntRange(0, 10000-1)*9;
			int sitechoice1 = rand_IntRange(0,2);
			int sitechoice2 = rand_IntRange(0,2);
			int vectorchoice = rand_IntRange(0,10000-1)*3;
			double distlist[9] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};

			for(int q = 0; q < 3 ; q++)
			{
				for(int qq = 0; qq < 3; qq++)
				{
					tempcoordinates1[q*3 + qq] = coordlist[molchoice1 + q*3 + qq] - coordlist[molchoice1 + sitechoice1*3 + qq] ;
					tempcoordinates2[q*3 + qq] = coordlist[molchoice2 + q*3 + qq] - coordlist[molchoice2 + sitechoice2*3 + qq] + vectorlist[vectorchoice + qq]*n*delr;

				}
			}
			float totalE = 0.0f;
			float solE = 0.0f;
			double HSenergy;
			for(int q = 0; q < 3; q++)
			{
				for(int qq = 0; qq <3; qq++)
				{
					double distancesqr = 0.0f;
					for(int qqq = 0; qqq < 3; qqq++)
					{
						distancesqr +=(tempcoordinates1[q*3 + qqq] - tempcoordinates2[qq*3 + qqq])*(tempcoordinates1[q*3 + qqq] - tempcoordinates2[qq*3 + qqq]);
					}
					double distance = sqrt(distancesqr);
					distlist[q*3 + qq] = distance;
					double disttemp = distance/delr;
					distindexlist[q*3+ qq] = round(disttemp);
					int rindex = ceil(disttemp);
					double rfactor = rindex - disttemp;
					//interpolation of solvation potential grid data
					solE += solE + solFactor*(wr[interactionlist[q*3 + qq]*ngrid + rindex - 1] + wr[interactionlist[q*3 + qq]*ngrid + rindex - 1]);
					if (distance >= sigmalist[interactionlist[q*3+qq]])
					{
						HSenergy = 0.0;
					}
					else
					{
						HSenergy = 1000000000.0;
					}
					totalE += HSenergy;
				}
			}
			for(int q = 0; q < 9; q++)
			{
				gr[interactionlist[q]*ngrid + distindexlist[q]] += exp(-(totalE + solE)/(kb*T))*n*n;
			}
//			cout << "mol1 choice = " << molchoice1 << endl;
//			cout << "mol2 choice = " << molchoice2 << endl;
//			cout << "site1 choice = " << sitechoice1 << endl;
//			cout << "site2 choice = " << sitechoice2 << endl;
//			cout << "vector " << endl;
//			for(int w = 0; w<3; w++)
//			{
//				cout << vectorlist[vectorchoice + w]<< endl;
//			}
//			cout << "chain 1"<< endl;
//			for(int w = 0; w < 9; w++)
//			{
//				cout << tempcoordinates1[w] << endl;
//			}
//			cout << "chain 2"<< endl;
//			for(int w = 0; w < 9; w++)
//			{
//				cout << tempcoordinates2[w] << endl;
//			}
//			cout << "distance list"<< endl;
//			for(int w = 0; w < 9; w++)
//			{
//				cout << distlist[w] << endl;
//			}
//			cout << "distance list"<< endl;
//			for(int w = 0; w < 9; w++)
//			{
//				cout << distindexlist[w] << endl;
//			}
		}
//		for(int w = 0;w<9;w++)
//		{
//			cout << "gr["<< interactionlist[w]*ngrid + distindexlist[w] << "] = "<< gr[interactionlist[w]*ngrid + distindexlist[w]] << endl;
//		}
	}
	cout << "code finished"<< endl;
	return 0;


}
