#ifndef UNTITLED_FUNCTIONS_H
#define UNTITLED_FUNCTIONS_H
#include<vector>
#include<iostream>
#include "GlobalData.h"
using namespace std;
#define SIZE 16
void displayArray(vector<vector<double>> c,int size);
void displayArray(vector<vector<double>> c);
vector<vector<double>> sumVectors(vector<vector<double>> a,vector<vector<double>> b);
struct LocalMatrixElem2{
    vector<vector<double>> H={
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<vector<double>> C={
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<double> P={0,0,0,0};
    LocalMatrixElem2(){};
};
vector<vector<double>> sumUpHglobal(vector<vector<double>> Hlocal,vector<vector<double>> HGlobal,int n,FEMGrid grid);//zalezne od siatki oraz typu elementu siatki!
vector<double> sumUpPglobal(vector<double> Plocal,int n, FEMGrid grid);
void displayVector(vector<double> arg);
double shapeFun(double a,double b,double signa,double signb);
void getCofactor(double A[SIZE][SIZE], double temp[SIZE][SIZE], int p, int q, int n);
int determinant(double A[SIZE][SIZE], int n);
void adjoint(double A[SIZE][SIZE],double adj[SIZE][SIZE]);
bool inverse(double A[SIZE][SIZE], double inverse[SIZE][SIZE]);
vector<vector<double>> trans(vector<vector<double>> A);
vector<vector<double>> inv(vector<vector<double>> A);
double det(vector<vector<double>> A);


void getCfactor(double M[SIZE][SIZE], double t[SIZE][SIZE], int p, int q, int n);
double DET(double M[SIZE][SIZE], int n);
void ADJ(double M[SIZE][SIZE],double adj[SIZE][SIZE]);
bool INV(double M[SIZE][SIZE], double inv[SIZE][SIZE]);
template<class T> void print(T A[SIZE][SIZE]);
double maxVal(vector<double> arg);
vector<vector<double>> gauss(vector<vector<double>> a);
#endif //UNTITLED_FUNCTIONS_H
