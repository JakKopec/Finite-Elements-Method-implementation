#ifndef UNTITLED_GLOBALDATA_H
#define UNTITLED_GLOBALDATA_H
#include "Node.h"
#include "Element.h"
#include "fstream"
#include "vector"
using namespace std;

class GlobalData
{
public:
    double W;//szeerokosc siatki
    double H;//wysokosc siatki
    int nH;//ilosc wezlow na H
    int nW;//ilosc wezlow na W
    int nE;//ilosc elementow
    int nN;//ilosc wezlow
    double dw; //W/nW-1
    double dh; //H/nH-1
    int heatConductionIndex;//wsp przewodzenia ciepla
    int schema;//schemat calkowania
    double denisity;
    double c;
    double alpha;
    double initialTemperature;
    double ambientTemperature;
    double simulationTime;
    double simulationStepTime;
    GlobalData();
};
class FEMGrid:public GlobalData
{
public:
    vector<Node> arrN;//nN
    vector<Element> arrE;//nE
    vector<vector<double>> CGlobal;
    vector<vector<double>> HGlobal;
    vector<double> PGlobal;
    vector<double> PFinal;
    vector<vector<double>> CdT;
    vector<double> CdTt0;
    vector<vector<double>> HFinal;
    vector<double> t0Vector;
    bool gauss(int n, vector<vector<double>> AB, vector<double> X);

    FEMGrid();
};


#endif //UNTITLED_GLOBALDATA_H
