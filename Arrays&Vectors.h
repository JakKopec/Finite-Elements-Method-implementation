#ifndef UNTITLED_ARRAYS_VECTORS_H
#define UNTITLED_ARRAYS_VECTORS_H
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
    LocalMatrixElem2(){};
};
vector<vector<double>> sumUpHglobal(vector<vector<double>> Hlocal,vector<vector<double>> HGlobal,int n,FEMGrid grid);//zalezne od siatki oraz typu elementu siatki!

#endif //UNTITLED_ARRAYS_VECTORS_H
