#include "GlobalData.h"
#include "Functions.h"

GlobalData::GlobalData() {
    fstream file;
    file.open("DATA.txt", std::ios_base::in | std::ios_base::out);
    vector<double> data;
    string stemp;
    double dtemp;
    while(!file.eof())
    {
        getline(file,stemp);
        dtemp=stod(stemp);
        data.push_back(dtemp);
    }
    W = data[0];
    H = data[1];
    nW = data[2];
    nH = data[3];
    heatConductionIndex = data[4];
    schema = data[5];
    nE = (nH - 1) * (nW - 1);
    nN = nH * nW;
    dw = W / (nW - 1);
    dh = H / (nH - 1);
    denisity = data[6];
    c = data[7];
    specificHeat=data[7];
    initialTemperature=data[8];
    ambientTemperature=data[9];
    alpha=data[10];
    simulationTime=data[11];
    simulationStepTime=data[12];
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
    cout << "\nGESTOSC:" << denisity;
    cout << "\nC:" << c;
    cout << "\nCIPELO WLASCIWE(?):" << specificHeat;
    cout << "\nTEMPERTAURA POCZATKOWA:" << initialTemperature;
    cout << "\nTEMPERATURA OTOCZENIA:" << ambientTemperature;
    cout << "\nALFA:" << alpha;
    cout << "\nCZAS SYMULACJI:" << simulationTime;
    cout << "\nDLUGOSC KROKU SYMULACJI:" << simulationStepTime << endl;
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
        arrN[i].t0 = this->initialTemperature;
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
        arrN[i].t0 = this->initialTemperature;
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


