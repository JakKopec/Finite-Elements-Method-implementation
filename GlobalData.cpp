/*STRUKTURA PLIKU
1- szerokosc
2- wysokosc
3- ilosc wezlow w poziomie
4- ilosc wezlow w pionie
5- wspolczynnik przewodzenia ciepla
6- schemat calkowania
7- ro
8- C
9- t0
10- alfa*/
#include "GlobalData.h"
#include "Functions.h"

GlobalData::GlobalData() {
    string data[12]; //jesli zostana dodane nowe dane to bedzie trzeba zaktualizowac rozmiar
    fstream file;
    file.open("DATA.txt", std::ios_base::in | std::ios_base::out);
    for (int a = 0; a < 12; a++) {
        getline(file, data[a]);
        //cout<<data[a]<<endl;
    }
    file.close();
    std::string::size_type sz;
    double temp = stod(data[0]);
    W = temp;
    temp = stod(data[1]);
    H = temp;
    nW = stoi(data[2]);
    nH = stoi(data[3]);
    heatConductionIndex = stoi(data[4]);
    schema = stoi(data[5]);
    nE = (nH - 1) * (nW - 1);
    nN = nH * nW;
    dw = W / (nW - 1);
    dh = H / (nH - 1);
    temp = stod(data[6]);
    ro = temp;
    temp = stod(data[7]);
    c=temp;
    temp = stod(data[8]);
    t0=temp;
    temp = stod(data[9]);
    alfa=temp;
    temp = stod(data[10]);
    talfa=temp;
    temp = stod(data[11]);
    dt=temp;
}

FEMGrid::FEMGrid() {
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
    cout << "\nT0:" << t0;
    cout << "\nALFA:" << alfa << endl;
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
    for (int i = 0; i < nN; i++) {
        arrN[i].t0 = this->t0;
        if (arrN[i].x == 0 || arrN[i].x == W || arrN[i].y == H || arrN[i].y == 0) {
            arrN[i].bc = true;
        }
        arrN[i].displayNode();
    }

    int temp = 0;
    int temp2 = 1;
    for (int i = 0; i < nE; i++) {
        arrE[i].nodes[0].operator=(arrN[i + temp]);
        arrE[i].nodes[1].operator=(arrN[i + nH + temp]);
        arrE[i].nodes[2].operator=(arrN[i + 1 + nH + temp]);
        arrE[i].nodes[3].operator=(arrN[i + 1 + temp]);
        arrE[i].displayElement();
        temp2++;
        if (temp2 % nH == 0) {
            temp++;
            temp2 = 1;
        }
    }
    cout << endl;

    for (int i = 0; i < nN; i++) {
        arrN[i].t0 = this->t0;
        if (arrN[i].x == 0 || arrN[i].x == W || arrN[i].y == H || arrN[i].y == 0) {
            arrN[i].bc = true;
        }
        arrN[i].displayNode();
    }

    CGlobal = vector<vector<double>>(nN, vector<double>(nN, 0));
    HGlobal = vector<vector<double>>(nN, vector<double>(nN, 0));
    for (int i = 0; i < this->nN; i++) {
        PGlobal.push_back(0);
    }
}


