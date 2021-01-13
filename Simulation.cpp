#include "Functions.h"
#include "ElemSolve.h"

void simulation(FEMGrid grid) {
    int iter = grid.simulationTime / grid.simulationStepTime;
    double temp;
    LocalMatrixElemData localMatrixElemData;
    //HFinal - H+(C/dtau), potem (H+(C/dtau))^-1
    //PFinal - P+((C/dtau)*t0)
    //tVector - wektor rozwiazan ukladu rownan
    //elem solve zwraca macierze lokalne
    for (int a = 0; a < grid.nE; a++) {
        //agregacja
        localMatrixElemData = elemSolve(grid.arrE[a], grid);
        grid.HGlobal = sumUpHglobal(localMatrixElemData.H, grid.HGlobal, a, grid);
        grid.CGlobal = sumUpHglobal(localMatrixElemData.C, grid.CGlobal, a, grid);
        grid.PGlobal = sumUpPglobal(localMatrixElemData.P, a, grid);
    }
    //displayArray(grid.HGlobal);
    //displayArray(grid.CGlobal);
    //displayVector(grid.PGlobal);

    vector<vector<double>> CdTau = vector<vector<double>>(grid.nN, vector<double>(grid.nN, 0));
    for (int i = 0; i < grid.nN; i++) {
        for (int j = 0; j < grid.nN; j++) {
            grid.HFinal[i][j] = grid.HGlobal[i][j] + (grid.CGlobal[i][j] / grid.simulationStepTime);
            CdTau[i][j] = grid.CGlobal[i][j] / grid.simulationStepTime;
        }
    }

    //POROWNANIE Z TEST CASE:

    /*cout << "Macierz H+Hbc:\n";
    displayArray(grid.HGlobal);
    cout << "Macierz C:\n";
    displayArray(grid.CGlobal);
    cout << "Wektor P:\n";
    displayVector(grid.PGlobal);

    cout << "Macierz H+C/dT:\n";
    displayArray(grid.HFinal);*/

    //________________________________

    grid.HFinal = gaussJordanEliminination(grid.HFinal, grid);


    for (int i = 0; i < iter; i++) {
        cout << "______________________Iteracja nr " << i + 1 << "______________________\n";
        //displayVector(grid.t0Vector);
        for (int i = 0; i < grid.nN; i++) {
            grid.PFinal[i]=grid.PGlobal[i];
            for (int j = 0; j < grid.nN; j++) {
                //grid.CdTt0[i] += (grid.t0Vector[j] * grid.CGlobal[i][j] / grid.simulationStepTime);
                grid.PFinal[i]+=CdTau[i][j]*grid.t0Vector[j];
            }
        }
        //cout << "Wektor P+C/dT*t0:\n";
        //displayVector(grid.PFinal);
        for (int i = 0; i < grid.nN;i++) {
            temp = 0;
            for (int j = 0; j < grid.nN; j++) {
                temp += grid.HFinal[i][j] * grid.PFinal[j];
            }
            grid.t0Vector[i] = temp;
        }

        cout << "Rozwiazanie ukladu rownan:\n";
        //displayVector(grid.t0Vector);
        cout << "Minimalna temperatura w zbiorze rozwiazan: ";
        cout.precision(6);
        cout.flush();
        cout << minVal(grid.t0Vector) << endl;
        cout << "Maksymalna temperatura w zbiorze rozwiazan: ";
        cout.precision(6);
        cout.flush();
        cout << maxVal(grid.t0Vector) << endl;

        grid.CdTt0 = vector<double>(grid.nN, 0);
    }
}