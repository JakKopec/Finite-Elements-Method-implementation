#ifndef UNTITLED_FUNCTIONS_H
#define UNTITLED_FUNCTIONS_H
#include<vector>
#include<iostream>
#include "GlobalData.h"
using namespace std;
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


#endif //UNTITLED_FUNCTIONS_H
