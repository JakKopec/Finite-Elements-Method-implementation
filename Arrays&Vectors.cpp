#include"Arrays&Vectors.h"
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
};
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

