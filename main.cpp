#include <iostream>
#include<vector>
#include<string>
#include<fstream>
#include <cmath>
#include<iomanip>
using namespace std;
/*STRUKTURA PLIKU
1- szerokosc
2- wysokosc
3- ilosc wezlow w poziomie
4- ilosc wezlow w pionie
5- wspolczynnik przewodzenia ciepla
6- schemat calkowania
7- ro
8- C
9- t0*/
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
    cout<<"WIERSZE:"<<c.size()<<endl<<"KOLUMNY:"<<c[0].size()<<endl;
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
class Node
{
public:
    double x;
    double y;
    int id;
    double t0;
    bool bc=false;
    void Multiply(double dw,double dh,int muldw,int muldh)//przypisanie wspolrzednych wezlom
    {
        this->x=(dw*muldw);
        this->y=(dh*muldh);
    }
    Node& operator=(const Node &a)
    {
        this->x=a.x;
        this->y=a.y;
        this->id=a.id;
        return *this;
    }
    void displayNode()
    {
        cout.precision(1);
        cout<<"Node ID:"<<id<<"\tx="<<this->x<<"\t\ty="<<this->y<<"\t\tt0="<<this->t0<<"\t\tflaga wb:"<<bc<<endl;
    }
    int returnId(){
        return id;
    }
    Node(){};
    Node(double a,double b)
    {
        this->x=a;
        this->y=b;
    }
    Node(double a,double b,int id)
    {
        this->x=a;
        this->y=b;
        this->id=id;
    }

};

class Element {
public:
    int elemID;
    static int staticElemID;
    vector<Node> nodes;
    vector<int> id;
    vector<vector<double>> H = {//zmiana ilosci punktow w elemencie - zmiana macierzy H i C
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> C = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };

    Element() {
        this->elemID = ++staticElemID;
        this->nodes.push_back(Node(0, 0));
        this->nodes.push_back(Node(0, 0));
        this->nodes.push_back(Node(0, 0));
        this->nodes.push_back(Node(0, 0));
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
    }

    Element(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
        this->elemID = ++staticElemID;
        this->nodes.push_back(Node(x1, y1));
        this->nodes.push_back(Node(x2, y2));
        this->nodes.push_back(Node(x3, y3));
        this->nodes.push_back(Node(x4, y4));
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
    }

    Element(Node a, Node b, Node c, Node d) {
        this->elemID = ++staticElemID;
        this->nodes.push_back(a);
        this->nodes.push_back(b);
        this->nodes.push_back(c);
        this->nodes.push_back(d);
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
    }

    void displayElement() {
        cout << "\nElement: " << this->elemID << "\tPunkty calkowania:\n";
        for (int i = 0; i < 4; i++) {
            cout.precision(1);
            cout << "Node id:" << this->nodes[i].id + 1 << "  x: " << this->nodes[i].x << "\t\ty: " << this->nodes[i].y
                 << endl;
        }
    }
};
int Element::staticElemID=0;
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

    GlobalData()
    {
        string data[9]; //jesli zostana dodane nowe dane to bedzie trzeba zaktualizowac rozmiar
        fstream file;
        file.open("DATA.txt", std::ios_base::in | std::ios_base::out);
        for(int a=0;a<9;a++)
        {
            getline(file, data[a]);
            //cout<<data[a]<<endl;
        }
        file.close();
        std::string::size_type sz;
        double temp=stod(data[0]);
        W=temp;
        double temp2=stod(data[1]);
        H=temp2;
        nW=stoi(data[2]);
        nH=stoi(data[3]);
        heatConductionIndex=stoi(data[4]);
        schema=stoi(data[5]);
        nE=(nH-1)*(nW-1);
        nN=nH*nW;
        dw=W/(nW-1);
        dh=H/(nH-1);
        double temp3=stod(data[6]);
        double temp4=stod(data[7]);
        double temp5=stod(data[8]);
        ro=temp3;
        c=temp4;
        t0=temp5;

    }
};
class FEMGrid:public GlobalData
{
public:
    vector<Node> arrN;//nN
    vector<Element> arrE;//nE
    vector<vector<double>> CGlobal;
    vector<vector<double>> HGlobal;
    FEMGrid() {
        GlobalData a;
        cout << "\tODCZYT Z PLIKU:";
        cout << "\nSZEROKOSC:" << W;
        cout << "\nWYSOKOSC:" << H;
        cout << "\nILOSC WEZLOW NA 1 ODCINKU W:" << nW;
        cout << "\nILOSC WEZLOW NA 1 ODCINKU H:" << nH;
        cout << "\nILOSC ELEMENTOW:" << nE;
        cout << "\ndw:" << dw;
        cout << "\ndh:" << dh;
        cout << "\nWSP PRZEW CIEPLA:" << heatConductionIndex;
        cout << "\nSCHEMAT CALKOWANIA:" << schema << "-PUNKTOWY";
        cout << "\nRO:" << ro;
        cout << "\nC:" << c;
        cout << "\nT0:" << t0 << endl;
        int index = 0;
        for (int a = 0; a < nN; a++) {
            arrN.push_back(Node(0, 0));
            arrE.push_back(Element());
            arrN[a].id = a;
        }
        int dHtimes = 0;
        int dWtimes = 0;
        int c = 0;
        for (int i = 0; i < nW; i++) {
            dHtimes = 0;
            for (int j = 0; j < nH; j++, c++) {
                arrN[c].Multiply(dw, dh, dWtimes, dHtimes);
                dHtimes++;
                //arrN[c].displayNode();
            }
            dWtimes++;
        }
        int temp = 0;
        int temp2= 1;
        for (int i = 0; i < nE; i++) {
            arrE[i].nodes[0].operator=(arrN[i + temp]);
            arrE[i].nodes[1].operator=(arrN[i + nH + temp]);
            arrE[i].nodes[2].operator=(arrN[i + 1 + nH + temp]);
            arrE[i].nodes[3].operator=(arrN[i + 1 + temp]);
            arrE[i].displayElement();
            temp2++;
            if(temp2%nH==0){temp++;temp2=1;}
        }
        cout<<endl;

        for(int i=0;i<nN;i++){
            arrN[i].t0=this->t0;
            if(arrN[i].x==0 || arrN[i].x==W || arrN[i].y==H || arrN[i].y==0){
                arrN[i].bc=true;
            }
            arrN[i].displayNode();
        }

        CGlobal=vector<vector<double>>(nN,vector<double>(nN,0));
        HGlobal=vector<vector<double>>(nN,vector<double>(nN,0));
        displayArray(CGlobal);
        displayArray(HGlobal);
    }


};
double twoVariablesFunction(double x,double y)
{
    //dla funkcji -5x^2y+2xy^2+10
    return (-5)*x*x*y+2*x*y*y+10;
}
double gaussQuadrature2D(int k){
    if(k==2 | k==3){
        double finalResult=0;
        double x,y=0;
        if(k==3){
            double points[3]={-0.774597,0,0.774597};
            double weights[3]={0.555555,0.888888,0.555555};
            for (int i=0;i<k;i++){
                for (int j=0;j<k;j++){
                    x=points[i];
                    y=points[j];
                    finalResult += twoVariablesFunction(x,y)*weights[i]*weights[j];
                }
            }
        }
        else if(k==2){
            double points[2]={-1/sqrt(3),1/sqrt(3)};
            double weights[2]={1,1};
            for (int i=0;i<k;i++){
                for (int j=0;j<k;j++){
                    x=points[i];
                    y=points[j];
                    finalResult += twoVariablesFunction(x,y)*weights[i]*weights[j];
                }
            }
        }

        return finalResult;
    }
    cout<<"Wrong argument!\n";
    return -1;
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
LocalMatrixElem2 elem2solve(Element b,FEMGrid grid) {
    double ksi[4] = {(-1 / sqrt(3)), (1 / sqrt(3)), (1 / sqrt(3)), (-1 / sqrt(3))};
    double eta[4] = {(-1 / sqrt(3)), (-1 / sqrt(3)), (1 / sqrt(3)), (1 / sqrt(3))};
    int weight[4] = {1, 1, 1, 1};
    double det = 0;
    vector<vector<double>> dNdKsi
            {
                    {(-0.25 * (1 - eta[0])), (0.25 * (1 - eta[0])), (0.25 * (1 + eta[0])), (-0.25 * (1 + eta[0]))},
                    {(-0.25 * (1 - eta[1])), (0.25 * (1 - eta[1])), (0.25 * (1 + eta[1])), (-0.25 * (1 + eta[1]))},
                    {(-0.25 * (1 - eta[2])), (0.25 * (1 - eta[2])), (0.25 * (1 + eta[2])), (-0.25 * (1 + eta[2]))},
                    {(-0.25 * (1 - eta[3])), (0.25 * (1 - eta[3])), (0.25 * (1 + eta[3])), (-0.25 * (1 + eta[3]))}
            };
    vector<vector<double>> dNdEta
            {
                    {(-0.25 * (1 - ksi[0])), (-0.25 * (1 + ksi[0])), (0.25 * (1 + ksi[0])), (0.25 * (1 - ksi[0]))},
                    {(-0.25 * (1 - ksi[1])), (-0.25 * (1 + ksi[1])), (0.25 * (1 + ksi[1])), (0.25 * (1 - ksi[1]))},
                    {(-0.25 * (1 - ksi[2])), (-0.25 * (1 + ksi[2])), (0.25 * (1 + ksi[2])), (0.25 * (1 - ksi[2]))},
                    {(-0.25 * (1 - ksi[3])), (-0.25 * (1 + ksi[3])), (0.25 * (1 + ksi[3])), (0.25 * (1 - ksi[3]))}
            };
    double jacobian[4] = {0, 0, 0, 0};
    double reversedJacobian[4] = {0, 0, 0, 0};
    vector<vector<double>> dNdX = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> dNdY = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> dNdXT = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> dNdYT = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> multipliedX = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> multipliedY = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> tempH = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> tempC = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> N = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> NT = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> C = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };


    for (int i = 0; i < 4; i++) {
        jacobian[0] += dNdKsi[0][i] * b.nodes[i].x;
        jacobian[1] += dNdEta[1][i] * b.nodes[i].x;
        jacobian[2] += dNdKsi[2][i] * b.nodes[i].y;
        jacobian[3] += dNdEta[3][i] * b.nodes[i].y;
    }
    det = (jacobian[0] * jacobian[3]) - (jacobian[1] * jacobian[2]);
    reversedJacobian[0] = jacobian[3] / det;
    reversedJacobian[1] = jacobian[1] / det;
    reversedJacobian[2] = jacobian[2] / det;
    reversedJacobian[3] = jacobian[0] / det;
    cout << "Element" << b.elemID << " \nDet:" << det << "\tJakobian:\n";
    for (int j = 0; j < 4; j++) {
        cout << jacobian[j] << "\t";
    }
    cout << endl;

    for (int a = 0; a < 4; a++) {
        dNdX[0][a] = (dNdKsi[0][a] * reversedJacobian[0] + dNdEta[0][a] * reversedJacobian[1]);//git
        dNdX[1][a] = (dNdKsi[1][a] * reversedJacobian[0] + dNdEta[1][a] * reversedJacobian[1]);//git
        dNdX[2][a] = (dNdKsi[2][a] * reversedJacobian[0] + dNdEta[2][a] * reversedJacobian[1]);//git
        dNdX[3][a] = (dNdKsi[3][a] * reversedJacobian[0] + dNdEta[3][a] * reversedJacobian[1]);//git
    }
    for (int a = 0; a < 4; a++) {
        dNdY[0][a] = (dNdKsi[0][a] * reversedJacobian[2] + dNdEta[0][a] * reversedJacobian[3]);//git
        dNdY[1][a] = (dNdKsi[1][a] * reversedJacobian[2] + dNdEta[1][a] * reversedJacobian[3]);//git
        dNdY[2][a] = (dNdKsi[2][a] * reversedJacobian[2] + dNdEta[2][a] * reversedJacobian[3]);//git
        dNdY[3][a] = (dNdKsi[3][a] * reversedJacobian[2] + dNdEta[3][a] * reversedJacobian[3]);//git
    }
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            dNdXT[a][b] = dNdX[b][a];
            dNdYT[a][b] = dNdY[b][a];
        }
    }


    for (int i = 0; i < 4; i++) {
        N[i][0] = (0.25 * (1 - ksi[i]) * (1 - eta[i]));
        N[i][1] = (0.25 * (1 + ksi[i]) * (1 - eta[i]));
        N[i][2] = (0.25 * (1 + ksi[i]) * (1 + eta[i]));
        N[i][3] = (0.25 * (1 - ksi[i]) * (1 + eta[i]));
    }

    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            NT[a][b] = N[b][a];
        }
    }
    for (int point = 1; point <= 4; point++) {

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                multipliedX[i][j] = dNdX[point - 1][j] * dNdXT[i][point - 1];
                multipliedY[i][j] = dNdY[point - 1][j] * dNdYT[i][point - 1];
            }
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                tempH[i][j] = grid.heatConductionIndex * det * (multipliedX[i][j] + multipliedY[i][j])* weight[point-1];
            }
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                tempC[i][j] = det * grid.ro * grid.c * (N[point - 1][j] * NT[i][point - 1]);
                cout<<tempC[i][j]<<"\t";
            }
            cout<<endl;
        }
        cout<<endl;
        b.H = sumVectors(b.H, tempH);
        b.C = sumVectors(b.C, tempC);

        /*int s;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                s = 0;
                for (int k = 0; k < 4; k++) { s += N[i][k] * NT[k][j]; }
                C[i][j] = (grid.ro * grid.c) * s;
                cout<<grid.ro<<"\t"<<grid.c<<"\t"<<s<<"\n";
            }*/
    }
    /*cout << endl;
    cout << "Element " << b.elemID << endl;
    cout << "dNdKsi:\n";
    displayArray(dNdKsi, 4);
    cout << "dNdEta:\n";
    displayArray(dNdEta, 4);
    cout << "dNdX:\n";
    displayArray(dNdX,4);
    cout << "dNdXT:\n";
    displayArray(dNdXT,4);
    cout << "dNdY:\n";
    displayArray(dNdY,4);
    cout << "dNdYT:\n";
    displayArray(dNdYT,4);
    cout << "multipliedX:\n";
    displayArray(multipliedX,4);
    cout << "multipliedY:\n";
    displayArray(multipliedY,4);*/
    cout << "Macierz H:\n";
    displayArray(b.H, 4);
    cout << "Macierz C:\n";
    displayArray(b.C, 4);


    LocalMatrixElem2 localMatrixElem2;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++) {
            localMatrixElem2.H[i][j] = b.H[i][j];
            localMatrixElem2.C[i][j] = b.C[i][j];
        }
    }
    return localMatrixElem2;
}





vector<vector<double>> elem3solve(Element b) {
    double sq = sqrt(15) / 5;
    double ksi[9] =
            {-sq, 0, sq, -sq, 0, sq, -sq, 0, sq};//do poprawy
    double eta[9] =
            {-sq, -sq, -sq, 0, 0, 0, sq, sq, sq};//do poprawy
    int weight[9] =
            {25 / 81, 40 / 81, 25 / 81, 40 / 81, 64 / 81, 40 / 81, 25 / 81, 40 / 81, 25 / 81};
    vector<vector<double>> jacobian{
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<vector<double>> reversedJacobian = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    double det[9] = {0};
    vector<vector<double>> dNdKsi{
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<vector<double>> dNdEta{
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};

    for (int i = 0; i < 9; i++) {
        dNdKsi[i][0] = -0.25 * (1 - eta[i]);
        dNdKsi[i][1] = 0.25 * (1 - eta[i]);
        dNdKsi[i][2] = 0.25 * (1 + eta[i]);
        dNdKsi[i][3] = -0.25 * (1 + eta[i]);

        dNdEta[i][0] = -0.25 * (1 - ksi[i]);
        dNdEta[i][1] = -0.25 * (1 - ksi[i]);
        dNdEta[i][2] = 0.25 * (1 + ksi[i]);
        dNdEta[i][3] = 0.25 * (1 + ksi[i]);
    }

    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 4; j++) {
            jacobian[i][0] += dNdKsi[0][j] * b.nodes[j].x;
            jacobian[i][1] += dNdEta[1][j] * b.nodes[j].x;
            jacobian[i][2] += dNdKsi[2][j] * b.nodes[j].y;
            jacobian[i][3] += dNdEta[3][j] * b.nodes[j].y;
        }
    }



    /*vector<double> dXdEta = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> dXdKsi = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> dYdEta = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> dYdKsi = {0, 0, 0, 0, 0, 0, 0, 0, 0};


    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 4; j++) {
            dXdKsi[i] += dNdKsi[i][j] * b.nodes[j].x;
            dYdKsi[i] += dNdKsi[i][j] * b.nodes[j].y;
            dXdEta[i] += dNdEta[i][j] * b.nodes[j].x;
            dYdEta[i] += dNdEta[i][j] * b.nodes[j].y;
        }
    }

    for (int i = 0; i < 9; i++) {
        jacobian[i][0] = dXdKsi[i];
        jacobian[i][1] = dYdKsi[i];
        jacobian[i][2] = dXdEta[i];
        jacobian[i][3] = dYdEta[i];
    }*/
    cout << "Jakobian\n";
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 4; j++) {
            cout << jacobian[i][j] << "\t\t";
        }
        cout << endl;
    }

    for (int i = 0; i < 9; i++) {
        det[i] = (jacobian[i][0] * jacobian[i][3]) - (jacobian[i][1] * jacobian[i][2]);
        cout<<"DET:"<<det<<endl;
        /*reversedJacobian[i][0] = (1 / det[i]) * jacobian[i][3];
        reversedJacobian[i][1] = (1 / det[i]) * -jacobian[i][1];
        reversedJacobian[i][2] = (1 / det[i]) * -jacobian[i][2];
        reversedJacobian[i][3] = (1 / det[i]) * jacobian[i][0];*/
    }



/*
        displayArray(reversedJacobian, 4);

        vector<vector<double>> dNdX{
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        vector<vector<double>> dNdY{
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        vector<vector<double>> H{
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},};

        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 4; j++) {
                dNdX[i][j] = (1 / det[i]) * (dNdKsi[i][j] * dYdEta[i] - dNdEta[i][j] * dYdKsi[i]);
                dNdY[i][j] = (1 / det[i]) * (-dXdEta[i] * dNdKsi[i][j] + dXdKsi[i] * dNdEta[i][j]);
            }
        }
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    b.H[i][j] += weight[i] * det[i] * 25 * (dNdX[i][j] * dNdX[i][k] + dNdY[i][j] * dNdY[i][k]);
                }
            }
        }
        cout << "\n\n\nMacierz H\n";
        displayArray(b.H, 4);*/
    return b.H;
}



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
        grid.HGlobal=sumUpHglobal(localMatrixElem2.H,grid.HGlobal,a,grid);
        grid.CGlobal=sumUpHglobal(localMatrixElem2.C,grid.CGlobal,a,grid);
    }
    displayArray(grid.HGlobal);
    displayArray(grid.CGlobal);

    //grid.arrE[0].H=elem3solve(grid.arrE[0]);
    //displayArray(grid.arrE[0].H,4);

}
