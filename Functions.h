#ifndef UNTITLED_FUNCTIONS_H
#define UNTITLED_FUNCTIONS_H
#include<vector>
#include<iostream>
#include "GlobalData.h"
using namespace std;
#define SIZE 16
void displayArray(vector<vector<double>> c,int size);
void displayArray(vector<vector<double>> c);
void displayArray(vector<vector<double>> c,int rows,int columns);
vector<vector<double>> sumVectors(vector<vector<double>> a,vector<vector<double>> b);
struct LocalMatrixElemData{
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
    LocalMatrixElemData(){};
};
vector<vector<double>> sumUpHglobal(vector<vector<double>> Hlocal,vector<vector<double>> HGlobal,int n,FEMGrid grid);//zalezne od siatki oraz typu elementu siatki!
vector<double> sumUpPglobal(vector<double> Plocal,int n, FEMGrid grid);
void displayVector(vector<double> arg);
double shapeFun(double a,double b,double signa,double signb);
double minVal(vector<double> arg);
double maxVal(vector<double> arg);
vector<vector<double>> gaussJordanEliminination(vector<vector<double>> a,FEMGrid femGrid);
void simulation(FEMGrid femGrid);
#endif //UNTITLED_FUNCTIONS_H
