#include <cmath>
#include"Functions.h"
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
void simulation(FEMGrid femGrid,int iteration)
{
    cout<<"Iteracja nr "<<iteration<<":\n";
}

#define SIZE 16
void getCfactor(double M[SIZE][SIZE], double t[SIZE][SIZE], int p, int q, int n) {
    int i = 0, j = 0;
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) //Copy only those elements which are not in given row r and column c: {
            if (r != p && c != q) {
                t[i][j++] = M[r][c]; //If row is filled increase r index and reset c index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
    }
}
double DET(double M[SIZE][SIZE], int n) {
    double D = 0;
    if (n == 1)
        return M[0][0];
    double t[SIZE][SIZE]; //store cofactors
    int s = 1;
    for (int f = 0; f < n; f++) {
//For Getting Cofactor of M[0][f] do getCfactor(M, t, 0, f, n); D += s * M[0][f] * DET(t, n - 1);
        s = -s;
    }
    return D;
}
void ADJ(double M[SIZE][SIZE],double adj[SIZE][SIZE]) {
    if (SIZE == 1) {
        adj[0][0] = 1;
        return;
    }
    int s = 1;
    double t[SIZE][SIZE];
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
//To get cofactor of M[i][j]
            getCfactor(M, t, i, j, SIZE);
            s = ((i + j) % 2 == 0) ? 1 : -1; //sign of adj[j][i] positive if sum of row and column indexes is even.
            adj[j][i] =
                    (s) * (DET(t, SIZE - 1)); //Interchange rows and columns to get the transpose of the cofactor matrix
        }
    }
}
bool INV(double M[SIZE][SIZE], double inv[SIZE][SIZE]) {
    int det = DET(M, SIZE);
    if (det == 0) {
        cout << "can't find its inverse";
        return false;
    }
    double adj[SIZE][SIZE];
    ADJ(M, adj);
    for (int i = 0; i < SIZE; i++) for (int j = 0; j < SIZE; j++) inv[i][j] = adj[i][j] / double(det);
    return true;
}
template<class T> void print(T A[SIZE][SIZE]) {
for (int i=0; i<SIZE; i++) { for (int j=0; j<SIZE; j++) cout << A[i][j] << " "; cout << endl; }
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
vector<vector<double>> gauss(vector<vector<double>> a) {
    cout.precision(3);
    vector<vector<double>> result= vector<vector<double>>(16, vector<double>(16, 0));
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
    cout << endl << "Input matrix:" << endl;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
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
    cout << endl;
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
    cout << endl;
    for (i = 0; i < n; ++i) {
        d = mat[i][i];
        for (j = 0; j < 2 * n; ++j) {
            mat[i][j] = mat[i][j] / d;
        }
    }
    for (i = 0; i < n; ++i) {
        for (j = n; j < 2 * n; ++j) {
            //cout << mat[i][j] << "\t";
            result[i][j-16]=mat[i][j];
        }
        //cout << endl;
    }
    displayArray(result);
    for (i = 0; i < n; ++i) {
        delete[] mat[i];
    }
    delete[] mat;
    return result;
}