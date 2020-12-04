#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "Arrays&Vectors.h"
#include "ElemSolve.h"
using namespace std;
int main() {
    FEMGrid grid;

    /*cout << "______________________________\n";
    double calka = gaussQuadrature2D(2);
    cout << endl << calka << endl << endl;
    cout << "______________________________\n";*/

    //element testowy
    /*vector<vector<double>> h={
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    Element test(0,4,4,0,0,0,6,6);
    h=elem3solve(test);*/

    //grid.arrE[0].H=elem2solve(grid.arrE[0]);
    //displayArray(grid.arrE[0].H,4);


    LocalMatrixElem2 localMatrixElem2;
    for (int a = 0; a < grid.nE; a++) {
        localMatrixElem2=elem2solve(grid.arrE[a],grid);
        //localMatrixElem2=elem3solve(grid.arrE[a],grid);
        grid.HGlobal=sumUpHglobal(localMatrixElem2.H,grid.HGlobal,a,grid);
        grid.CGlobal=sumUpHglobal(localMatrixElem2.C,grid.CGlobal,a,grid);
    }
    displayArray(grid.HGlobal);
    cout<<endl<<endl;
    displayArray(grid.CGlobal);

    //grid.arrE[0].H=elem3solve(grid.arrE[0]);
    //displayArray(grid.arrE[0].H,4);

}
