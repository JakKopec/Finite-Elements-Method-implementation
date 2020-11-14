#include <iostream>
#include<vector>
#include<string>
#include<fstream>
#include <cmath>
#include <iomanip>
using namespace std;
/*STRUKTURA PLIKU
1- szerokosc
2- wysokosc
3- ilosc wezlow w poziomie
4- ilosc wezlow w pionie*/
class Elem4;
class Element;
double elem4solve(Elem4 a,Element b);
class Node
{
public:
    double x;
    double y;
    void Multiply(double dw,double dh,double muldw,double muldh)//przypisanie wspolrzednych wezlom
    {
        this->x=(dw*muldw);
        this->y=(dh*muldh);
    }
    void overwrite(Node a)
    {
        this->x=a.x;
        this->y=a.y;
    }
    Node& operator=(const Node &a)
    {
        this->x=a.x;
        this->y=a.y;
        return *this;
    }
    void displayNode()
    {
        cout.precision(3);
        cout<<"x="<<this->x<<"\t\ty="<<this->y<<endl;
    }
    Node(){};
    Node(double a,double b)
    {
        this->x=a;
        this->y=b;
    }
};
static int staticelemID=0;
class Element
{
    friend double elem4solve(Elem4 a,Element b);
public:
    int elemID;
    vector <Node> nodes;
    vector <int> id;
    Element()
    {
        this->elemID=++staticelemID;
        this->nodes.push_back(Node(0,0));
        this->nodes.push_back(Node(0,0));
        this->nodes.push_back(Node(0,0));
        this->nodes.push_back(Node(0,0));
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
    }
    Element(double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4)
    {
        this->elemID=++staticelemID;
        this->nodes.push_back(Node(x1,y1));
        this->nodes.push_back(Node(x2,y2));
        this->nodes.push_back(Node(x3,y3));
        this->nodes.push_back(Node(x4,y4));
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
        this->id.push_back(0);
    }
    void displayElement()
    {
        for(int i=0;i<4;i++)
        {
            cout<<this->nodes[i].x<<"\t"<<this->nodes[i].y<<endl;
        }
        return;
    }

    void overwriteID(int a)//przy tworzeniu siatki konieczna jest kontrola polozenia elementu(wiersz/kolumna)
    {
        this->id[0]=this->elemID;
        this->id[1]=this->elemID+a;
        this->id[2]=this->elemID+a+1;
        this->id[3]=this->elemID+1;
    }
};
class GlobalData
{
public:
    int W;//szeerokosc siatki
    int H;//wysokosc siatki
    int nH;//ilosc wezlow na H
    int nW;//ilosc wezlow na W
    int nE;//ilosc elementow
    int nN;//ilosc wezlow
    double dw; //W/nW-1
    double dh; //H/nH-1

    GlobalData()
    {
        string data[4]; //jesli zostana dodane nowe dane to bedzie trzeba zaktualizowac rozmiar
        fstream file;
        file.open("DATA.txt", std::ios_base::in | std::ios_base::out);
        for(int a=0;a<4;a++)
        {
            getline(file, data[a]);
            //cout<<data[a]<<endl;
        }
        file.close();
        std::string::size_type sz;
        W=stoi(data[0]);
        H=stoi(data[1]);
        nW=stoi(data[2]);
        nH=stoi(data[3]);
        nE=(nH-1)*(nW-1);
        nN=nH*nW;
        dw=(double)W/(nW-1);
        dh=(double)H/(nH-1);
        cout<<"\tODCZYT Z PLIKU:";
        cout<<"\nSZEROKOSC:"<<W;
        cout<<"\nWYSOKOSC:"<<H;
        cout<<"\nILOSC WEZLOW NA 1 ODCINKU W:"<<nW;
        cout<<"\nILOSC WEZLOW NA 1 ODCINKU H:"<<nH;
        cout<<"\nILOSC ELEMENTOW:"<<nE;
        cout<<"\ndw:"<<dw;
        cout<<"\ndh:"<<dh<<endl<<endl;
    }
};

class FEMGrid:public GlobalData
{
public:
    vector<Node> arrN;//nN
    vector<Element> arrE;//nE
    friend class Elem4;

    FEMGrid() {
        GlobalData a;
        int index = 0;
        for (int a = 0; a < nN; a++) {
            arrN.push_back(Node(0, 0));
            arrE.push_back(Element());
        }

        int temp = 0;//sposob na "oszukanie" multiply aby wezly.x z pierwszej kolumny byly wyzerowane
        for (int a = 0, b = 0; a < this->nN; a++, b++)//Przypisanie wspolrzednych wezlom
        {
            if (a / nH == temp) {
                //cout<<"MUL: temp="<<temp<<"\tb="<<b<<"\n"; //testy
                arrN[a].Multiply(this->dw, this->dh, temp, b);
            } else if (a / nH != temp) {
                temp++;
                b = 0;
                //cout<<"MUL: temp="<<temp<<"\tb="<<b<<"\n"; //testy
                arrN[a].Multiply(this->dw, this->dh, temp, b);
            }
            cout << "\nNode " << a + 1 << "\t\t";
            arrN[a].displayNode();//wypisanie wezlow
        }
        //Przypisanie elementom siatki poszczegolne wezly
        int numElements = 0;
        int nextColumn = 0;//zmienna pomagajaca w "przeskakiwaniu" do nowej kolumny elementow siatki
        for (int b = 0; b < this->nW - 1; b++,nextColumn++) {
            //cout<<"\nNEXT\n"<<nextColumn;
            for (int c = 0; c < this->nH - 1; c++, numElements++) {
                /*cout<<"Element\t"<<numElements
                <<"\tWezly:\t"<<numElements+nextColumn+1<<"\t"
                <<numElements + this->nH + nextColumn+1<<"\t"
                <<numElements + this->nH + 1 + nextColumn+1<<"\t"
                <<numElements + 1 + nextColumn+1<<endl;*///testy
                //wpisanie wspolrzednych wezlow do elementow

                arrE[numElements].nodes[0].operator=(arrN[numElements+nextColumn]);
                arrE[numElements].nodes[1].operator=(arrN[numElements+nextColumn+nH]);
                arrE[numElements].nodes[2].operator=(arrN[numElements+nextColumn+nH+1]);
                arrE[numElements].nodes[3].operator=(arrN[numElements+nextColumn+1]);
                displayElementsData(numElements);
            }
        }
    }
    void displayElementsData(int index)
    {
        cout<<"\nElement "<<arrE[index].elemID<<endl;
        arrE[index].nodes[0].displayNode();
        arrE[index].nodes[1].displayNode();
        arrE[index].nodes[2].displayNode();
        arrE[index].nodes[3].displayNode();
    }
};


double twoVariablesFunction(double x,double y)
{
    //dla funkcji -5x^2y+2xy^2+10
    return (-5)*x*x*y+2*x*y*y+10;
}
double gaussQuadrature2D(const int k){

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
void displayArray(vector<vector<double>> c) {
    cout.precision(3);
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            cout << c[a][b] << "\t";
        }
        cout << endl;
    }
    cout<<endl<<endl;
}

class Elem4 {
public:
    double ksi[4] = {(-1 / sqrt(3)), (1 / sqrt(3)), (1 / sqrt(3)), (-1 / sqrt(3))};
    double eta[4] = {(-1 / sqrt(3)), (-1 / sqrt(3)), (1 / sqrt(3)), (1 / sqrt(3))};
    int weight[4] = {1, 1, 1, 1};
    double jacobian[2][2] = {0};
    double reversedJacobian[2][2] = {0};
    double det = 0;
    vector<vector<double>> H = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
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

    /*vector<double> dNdKsi
            {
                    {(-0.25*(1-eta[0])),(0.25*(1-eta[0])),(0.25*(1+eta[0])),(-0.25*(1+eta[0]))}
            };
    vector<double> dNdEta
            {
                    {(-0.25*(1-ksi[0])),(-0.25*(1+ksi[0])),(0.25*(1+ksi[0])),(0.25*(1-ksi[0]))}
            };*/
    vector<vector<double>> elem4solve(int point, Element b, vector<double> x, vector<double> y) {
        /*cout<<"ETA:\n";
        for(int a=0;a<4;a++) {
            for (int b = 0; b < 4; b++) {
                cout << dNdEta[a][b] << "\t";
            }
            cout << endl;
        }
        cout<<"KSI:\n";
        for(int c=0;c<4;c++) {
            for (int d = 0; d < 4; d++) {
                cout << dNdKsi[c][d] << "\t";
            }
            cout << endl;
        }*/

        /*
        cout<<"dX/dEta\t\tdX/dKsi\n";
        for(int a=0;a<4;a++)
        {
            cout<<dNdEta[a]<<"\t"<<dNdKsi[a]<<endl;
        }*/
        double jacobian[4] = {0, 0, 0, 0};
        double reversedJacobian[4] = {0, 0, 0, 0};
        double det = 0;
        double i11 = 0;
        double i12 = 0;
        double i21 = 0;
        double i22 = 0;
        for (int i = 0; i < 4; i++) {
            i11 += dNdKsi[point - 1][i] * b.nodes[i].x;
            i21 += dNdEta[point - 1][i] * b.nodes[i].x;
            i12 += dNdKsi[point - 1][i] * b.nodes[i].y;
            i22 += dNdEta[point - 1][i] * b.nodes[i].y;
            /*i11+=dNdKsi[i]*b.nodes[i].x;
            i21+=dNdEta[i]*b.nodes[i].x;
            i12+=dNdKsi[i]*b.nodes[i].y;
            i22+=dNdEta[i]*b.nodes[i].y;*/
        }
        jacobian[0] = i11;
        jacobian[1] = i12;
        jacobian[2] = i21;
        jacobian[3] = i22;
        /*cout<<"JACOBIAN ARRAY:\n";
        for(int j=0;j<4;j++){
           cout<<jacobian[j]<<"\t";
        }
        cout<<endl;*/
        det = (i11 * i22) - (i12 * i21);
        reversedJacobian[0] = jacobian[3] / det;
        reversedJacobian[1] = jacobian[1] / det;
        reversedJacobian[2] = jacobian[2] / det;
        reversedJacobian[3] = jacobian[0] / det;
        /*cout << "dNdKsi:\n";
        displayArray(dNdKsi);
        cout << "dNdEta:\n";
        displayArray(dNdEta);*/

        vector<vector<double>> dNdX = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
        for (int a = 0; a < 4; a++) {
            dNdX[0][a] = (dNdKsi[0][a] * reversedJacobian[0] + dNdEta[0][a] * reversedJacobian[1]);//git
            dNdX[1][a] = (dNdKsi[1][a] * reversedJacobian[0] + dNdEta[1][a] * reversedJacobian[1]);//git
            dNdX[2][a] = (dNdKsi[2][a] * reversedJacobian[0] + dNdEta[2][a] * reversedJacobian[1]);//git
            dNdX[3][a] = (dNdKsi[3][a] * reversedJacobian[0] + dNdEta[3][a] * reversedJacobian[1]);//git
        }

        vector<vector<double>> dNdY = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
        for (int a = 0; a < 4; a++) {
            dNdY[0][a] = (dNdKsi[0][a] * reversedJacobian[2] + dNdEta[0][a] * reversedJacobian[3]);//git
            dNdY[1][a] = (dNdKsi[1][a] * reversedJacobian[2] + dNdEta[1][a] * reversedJacobian[3]);//git
            dNdY[2][a] = (dNdKsi[2][a] * reversedJacobian[2] + dNdEta[2][a] * reversedJacobian[3]);//git
            dNdY[3][a] = (dNdKsi[3][a] * reversedJacobian[2] + dNdEta[3][a] * reversedJacobian[3]);//git
        }

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
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                dNdXT[a][b] = dNdX[b][a];
                dNdYT[a][b] = dNdY[b][a];
            }
        }
        /*cout << "dNdX:\n";
        displayArray(dNdX);
        cout << "dNdXT:\n";
        displayArray(dNdXT);
        cout << "dNdY:\n";
        displayArray(dNdY);
        cout << "dNdYT:\n";
        displayArray(dNdYT);*/

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


        //WERSJA NR 2
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                multipliedX[i][j] = dNdX[point - 1][j] * dNdXT[i][point - 1];
            }
        }
        //displayArray(multipliedX);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                multipliedY[i][j] = dNdY[point - 1][j] * dNdYT[i][point - 1];
            }
        }
        //displayArray(multipliedY);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                H[i][j] = 30 * det * (multipliedX[i][j] + multipliedY[i][j]);//*det*30;
            }
        }

        cout << point << " punkt calkowania - macierz H:\n";
        displayArray(H);


        return H;
    }
};
vector<vector<double>> sumVectors4x4(vector<vector<double>> a,vector<vector<double>> b,vector<vector<double>> result){
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i][j]=a[i][j]+b[i][j];
        }
    }
    return result;
};

int main()
{
    //FEMGrid a;

    /*cout<<"______________________________\n";
    double calka=gaussQuadrature2D(2);
    cout<<endl<<calka<<endl<<endl;
    cout<<"______________________________\n";*/
    Elem4 b;
    Element c(0,4,4,0,0,0,6,6);
    vector<double> x={0,4,4,0};
    vector<double> y={0,0,6,6};

    vector<vector<double>> result={
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<vector<double>> temp{
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };

    for(int z=1;z<=4;z++) {
        temp = b.elem4solve(z, c, x, y);
        result = sumVectors4x4(temp, result, result);
        temp = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
    }

    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            cout << result[a][b] << "\t";
        }
        cout << endl;
    }
    cout<<endl<<endl;



    return 0;
}

