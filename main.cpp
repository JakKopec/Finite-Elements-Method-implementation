#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "Functions.h"
#include "Elem2Solve.h"
#include "matrix.h"

using namespace std;
int main() {
    FEMGrid grid;
    simulation(grid);
    /*LocalMatrixElem2 localMatrixElem2;
    //Element test(0,4,4,0,0,0,6,6);
    //localMatrixElem2=elem2solve(test,grid);
    for (int a = 0; a < grid.nE; a++) {
        localMatrixElem2 = elem2solve(grid.arrE[a], grid);
        //localMatrixElem2=elem3solve(grid.arrE[a],grid);
        grid.HGlobal = sumUpHglobal(localMatrixElem2.H, grid.HGlobal, a, grid);
        grid.CGlobal = sumUpHglobal(localMatrixElem2.C, grid.CGlobal, a, grid);
        grid.PGlobal = sumUpPglobal(localMatrixElem2.P, a, grid);
    }
    cout << "Macierz H+Hbc:\n";
    displayArray(grid.HGlobal);
    cout << "Macierz C:\n";
    displayArray(grid.CGlobal);
    cout << "Wektor P:\n";
    displayVector(grid.PGlobal);

    for (int i = 0; i < grid.nN; i++) {
        for (int j = 0; j < grid.nN; j++) {
            grid.HFinal[i][j] = grid.HGlobal[i][j] + grid.CGlobal[i][j] / grid.simulationStepTime;
        }
    }
    for (int i = 0; i < grid.nN; i++) {
        for (int j = 0; j < grid.nN; j++) {
            grid.CdTt0[i] += (grid.t0Vector[j] * grid.CGlobal[i][j] / grid.simulationStepTime);
        }
        grid.PFinal[i] = grid.PGlobal[i] + grid.CdTt0[i];
    }
    cout << "Macierz H+C/dT:\n";
    displayArray(grid.HFinal);
    cout << "Wektor P+C/dT*t0:\n";
    displayVector(grid.PFinal);
    cout << endl;

    vector<vector<double>> inverted = vector<vector<double>>(16, vector<double>(16, 0));
    vector<double> tvector = vector<double>(grid.nN,0);


    inverted=gaussJordanEliminination(grid.HFinal);
    double temp=0;
    for(int a=0;a<grid.nN;a++){
        temp=0;
        for(int b=0;b<grid.nN;b++){
            temp+=inverted[a][b]*grid.PFinal[b];
        }
        tvector[a]=temp;
    }
    displayVector(tvector);
    temp=maxVal(tvector);
    cout<<"Maksymalna temperatura w zbiorze rozwiazan: "<<temp<<endl;*/
}

         /*vector<vector<double>> aaa = {
            {1,18,9,9,5,3},
            {0,1,4,6,2,8},
            {48,9,4,7,5,1},
            {5,6,5,2,1,5},
            {8,9,77,8,5,1},
            {21,5,84,8,48,7}
    };*/