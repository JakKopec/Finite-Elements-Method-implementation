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
    double ro;
    double c;
    double t0; //temperatura poczatkowa w wezle
    double alfa;
    GlobalData();
};
class FEMGrid:public GlobalData
{
public:
    vector<Node> arrN;//nN
    vector<Element> arrE;//nE
    vector<vector<double>> CGlobal;
    vector<vector<double>> HGlobal;
    FEMGrid();
};


#endif //UNTITLED_GLOBALDATA_H
