#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "Functions.h"
#include "ElemSolve.h"
using namespace std;
int main() {
    FEMGrid grid;
    LocalMatrixElem2 localMatrixElem2;
    for (int a = 0; a < grid.nE; a++) {
        localMatrixElem2=elem2solve(grid.arrE[a],grid);
        //localMatrixElem2=elem3solve(grid.arrE[a],grid);
        grid.HGlobal=sumUpHglobal(localMatrixElem2.H,grid.HGlobal,a,grid);
        grid.CGlobal=sumUpHglobal(localMatrixElem2.C,grid.CGlobal,a,grid);
    }
    cout<<"Macierz H:\n";
    displayArray(grid.HGlobal);
    cout<<endl<<endl;
    cout<<"Macierz C:\n";
    displayArray(grid.CGlobal);

    //grid.arrE[0].H=elem3solve(grid.arrE[0]);
    //displayArray(grid.arrE[0].H,4);
}
