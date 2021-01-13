#include <cmath>
#include <iomanip>
#include"Functions.h"
#include "ElemSolve.h"

void displayArray(vector<vector<double>> c,int size) {
    cout.precision(3);
    for (int a = 0; a < size; a++) {
        for (int b = 0; b < size; b++) {
            cout << c[a][b] << "\t";
        }
        cout << endl;
    }
    cout<<endl;
}
void displayArray(vector<vector<double>> c) {
    cout.precision(3);
    //cout<<"WIERSZE:"<<c.size()<<endl<<"KOLUMNY:"<<c[0].size()<<endl;
    for (int a = 0; a < c.size(); a++) {
        for (int b = 0; b < c[0].size(); b++) {
            cout << c[a][b] << "\t";
        }
        cout << endl;
    }
    cout<<endl;
}
void displayArray(vector<vector<double>> c,int rows,int columns){
    cout.precision(3);
    for (int a = 0; a < rows; a++) {
        for (int b = 0; b < columns; b++) {
            cout << c[a][b] << "\t";
        }
        cout << endl;
    }
    cout<<endl;
}

vector<vector<double>> sumVectors(vector<vector<double>> a,vector<vector<double>> b){
    vector<vector<double>> result={
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a.size(); j++) {
            result[i][j]=a[i][j]+b[i][j];
        }
    }
    return result;
}
vector<vector<double>> sumUpHglobal(vector<vector<double>> Hlocal,vector<vector<double>> HGlobal,int n,FEMGrid grid)//zalezne od siatki oraz typu elementu siatki!
{
    int ind1, ind2 = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ind1 = grid.arrE[n].nodes[i].returnId();//cout<<ind1<<endl;
            ind2 = grid.arrE[n].nodes[j].returnId();//cout<<ind2<<endl;
            HGlobal[ind1][ind2] += Hlocal[i][j];
        }
    }
    //displayArray(HGlobal);
    return HGlobal;
}
vector<double> sumUpPglobal(vector<double> Plocal,int n, FEMGrid grid) {
    int ind;
    for (int i = 0; i < 4; i++) {
        ind=grid.arrE[n].nodes[i].returnId();
        //cout<<ind<<"\t";
        grid.PGlobal[ind]+=Plocal[i];
    }
    //cout<<endl;
    return grid.PGlobal;
}
double shapeFun(double a,double b,double signa,double signb){
    return 0.25*(1+signa*a)*(1+signb*b);
}
void displayVector(vector<double> arg){
    cout.precision(5);
    for(int a=0;a<arg.size();a++){
        cout<<arg[a]<<"\t";
    }
    cout<<endl<<endl;
}
double maxVal(vector<double> arg)
{
    double temp=arg[0];
    for(int a=1;a<arg.size();a++){
        if(arg[a]>temp){
            temp=arg[a];
        }
    }
    return temp;
}
double minVal(vector<double> arg)
{
    double temp=arg[0];
    for(int a=1;a<arg.size();a++){
        if(arg[a]<temp){
            temp=arg[a];
        }
    }
    return temp;
}
vector<vector<double>> gaussJordanEliminination(vector<vector<double>> a,FEMGrid femGrid) {
    cout.precision(3);
    vector<vector<double>> result= vector<vector<double>>(femGrid.nN, vector<double>(femGrid.nN, 0));
    int i = 0, j = 0, k = 0, n = 0;
    float **mat = NULL;
    float d = 0.0;
    n=a.size();
    mat = new float *[2 * n];
    for (i = 0; i < 2 * n; ++i) {
        mat[i] = new float[2 * n]();
    }
    for (i = 0; i < n; ++i) {
        for (j = 0;j < n;++j) {
            mat[i][j]=a[i][j];
        }
    }
    for (i = 0; i < n; ++i) {
        for (j = 0; j < 2 * n; ++j) {
            if (j == (i + n)) {
                mat[i][j] = 1;
            }
        }
    }
    for (i = n; i > 1; --i) {
        if (mat[i - 1][1] < mat[i][1]) {
            for (j = 0; j < 2 * n; ++j) {
                d = mat[i][j];
                mat[i][j] = mat[i - 1][j];
                mat[i - 1][j] = d;
            }
        }
    }
    for (i = 0; i < n; ++i) {
        for (j = 0; j < 2 * n; ++j) {
            if (j != i) {
                d = mat[j][i] / mat[i][i];
                for (k = 0; k < n * 2; ++k) {
                    mat[j][k] -= mat[i][k] * d;
                }
            }
        }
    }
    for (i = 0; i < n; ++i) {
        d = mat[i][i];
        for (j = 0; j < 2 * n; ++j) {
            mat[i][j] = mat[i][j] / d;
        }
    }
    for (i = 0; i < n; ++i) {
        for (j = n; j < 2 * n; ++j) {
            //cout << mat[i][j] << "\t";
            result[i][j-result.size()]=mat[i][j];
        }
        //cout << endl;
    }
    //displayArray(result);
    for (i = 0; i < n; ++i) {
        delete[] mat[i];
    }
    delete[] mat;
    return result;
}
void simulation(FEMGrid grid)
{
    int iter=grid.simulationTime/grid.simulationStepTime;
    LocalMatrixElemData localMatrixElemData;
    vector<vector<double>> invertedH;
    vector<double> tVector = vector<double>(grid.nN,0);
    for(int i=0;i<iter;i++) {
        cout << "______________________Iteracja nr " << i << "______________________\n";
        for (int a = 0; a < grid.nE; a++) {
            localMatrixElemData = elemSolve(grid.arrE[a], grid);
            grid.HGlobal = sumUpHglobal(localMatrixElemData.H, grid.HGlobal, a, grid);
            grid.CGlobal = sumUpHglobal(localMatrixElemData.C, grid.CGlobal, a, grid);
            grid.PGlobal = sumUpPglobal(localMatrixElemData.P, a, grid);
        }

        /*cout << "Macierz H+Hbc:\n";
        displayArray(grid.HGlobal);
        cout << "Macierz C:\n";
        displayArray(grid.CGlobal);
        cout << "Wektor P:\n";
        displayVector(grid.PGlobal);*/

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
        //cout<<grid.PFinal[0]<<endl;
        /*cout << "Macierz H+C/dT:\n";
        displayArray(grid.HFinal);
        cout << "Wektor P+C/dT*t0:\n";
        displayVector(grid.PFinal);*/
        invertedH=vector<vector<double>>(grid.nN, vector<double>(grid.nN, 0));
        invertedH=gaussJordanEliminination(grid.HFinal,grid);
        //cout << "Odwrocone H\n";
        //displayArray(invertedH);
        double temp=0;
        for(int a=0;a<grid.nN;a++){
            temp=0;
            for(int b=0;b<grid.nN;b++){
                temp+=invertedH[a][b]*grid.PFinal[b];
            }
            tVector[a]=temp;
        }
        cout << "Rozwiazanie ukladu rownan:\n";
        //displayVector(tVector);
        cout<<"Minimalna temperatura w zbiorze rozwiazan: ";cout.precision(6);cout.flush();cout<<minVal(tVector)<<endl;
        cout<<"Maksymalna temperatura w zbiorze rozwiazan: ";cout.precision(6);cout.flush();cout<<maxVal(tVector)<<endl;

        invertedH.clear();
        tVector.clear();
        grid.HGlobal.clear();
        grid.CGlobal.clear();
        grid.PGlobal.clear();
        grid.HFinal.clear();
        grid.PFinal.clear();
        grid.CdTt0.clear();
        grid.HGlobal=vector<vector<double>>(grid.nN, vector<double>(grid.nN, 0));
        grid.CGlobal=vector<vector<double>>(grid.nN, vector<double>(grid.nN, 0));
        grid.PGlobal=vector<double>(grid.nN,0);
        grid.HFinal=vector<vector<double>>(grid.nN, vector<double>(grid.nN, 0));
        grid.PFinal=vector<double>(grid.nN,0);
        grid.CdTt0=vector<double>(grid.nN,0);
        for(int b=0;b<grid.nN;b++)
        {
            grid.t0Vector[b]=tVector[b];
        }
        //invertedH = vector<vector<double>>(grid.nN, vector<double>(grid.nN, 0));
        tVector = vector<double>(grid.nN,0);
    }
}