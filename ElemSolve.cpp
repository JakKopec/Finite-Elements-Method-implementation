#include "ElemSolve.h"
LocalMatrixElemData elemSolve(Element b, FEMGrid grid) {
    double schema2 = grid.schema * grid.schema;
    vector<double> ksi(schema2, 0);
    vector<double> eta(schema2, 0);
    vector<double> weights(schema2, 0);
    vector<double> det(schema2, 0);
    vector<vector<double>> dNdKsi(schema2, vector<double>(4, 0));
    vector<vector<double>> dNdEta(schema2, vector<double>(4, 0));

    double sq,sq1,sq2,w1,w2;
    if(grid.schema==2) {
        /*ksi = {(-1 / sqrt(3)), (1 / sqrt(3)), (1 / sqrt(3)), (-1 / sqrt(3))};
        eta = {(-1 / sqrt(3)), (-1 / sqrt(3)), (1 / sqrt(3)), (1 / sqrt(3))};
        for(int i=0;i<schema2;i++){
            weights.push_back(1);
        }*/
        return elem2solve(b, grid);
    }
    else if (grid.schema==3) {
        sq = 0.2 * sqrt(15);
        w1 = 5.0 / 9;
        w2 = 8.0 / 9;
        weights = {w1 * w1, w1 * w2, w1 * w1, w1 * w2, w2 * w2, w2 * w1, w1 * w1, w1 * w2, w1 * w1};
        ksi = {-sq, 0, sq, -sq, 0, sq, -sq, 0, sq};
        eta = {-sq, -sq, -sq, 0, 0, 0, sq, sq, sq};
    }
    else if (grid.schema==4) {
        w1 = 1.0 / 36 * (18 - sqrt(30));
        w2 = 1.0 / 36 * (18 + sqrt(30));
        weights = {w1 * w1, w1 * w2, w1 * w2, w1 * w1,
                   w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                   w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                   w1 * w1, w1 * w2, w1 * w2, w1 * w1};
        sq1 = 1.0 / 35 * sqrt(525 + 70 * sqrt(30));
        sq2 = 1.0 / 35 * sqrt(525 - 70 * sqrt(30));
        for (int i = 0; i < 16; i++) {
            ksi.push_back(0);
            eta.push_back(0);
        }
        for (int i = 0; i < 4; i++) {
            ksi[i * 4] = (-sq1);
            ksi[i * 4 + 1] = (-sq2);
            ksi[i * 4 + 2] = (sq2);
            ksi[i * 4 + 3] = (sq1);
            eta[i] = (-sq1);
            eta[i + 4] = (-sq2);
            eta[i + 8] = sq2;
            eta[i + 12] = sq1;
        }
    }
    vector<vector<double>> jacobian(schema2, vector<double>(schema2, 0));
    vector<vector<double>> reversedJacobian(schema2, vector<double>(schema2, 0));
    vector<vector<double>> dNdX(schema2, vector<double>(4, 0));
    vector<vector<double>> dNdY(schema2, vector<double>(4, 0));
    vector<vector<double>> dNdXT(4, vector<double>(schema2, 0));
    vector<vector<double>> dNdYT(4, vector<double>(schema2, 0));
    //vector<vector<double>> multipliedX(schema2, vector<double>(schema2, 0));
    //vector<vector<double>> multipliedY(schema2, vector<double>(schema2, 0));
    vector<vector<double>> tempH(4, vector<double>(4, 0));
    vector<vector<double>> tempC(4, vector<double>(4, 0));
    vector<vector<double>> N(schema2, vector<double>(4, 0));
    vector<vector<double>> NT(4, vector<double>(schema2, 0));
    vector<vector<double>> C(4, vector<double>(4, 0));
    vector<vector<double>> HBc(4, vector<double>(4, 0));
    vector<double> NBc1 = {0, 0, 0, 0};
    vector<double> NBc2 = {0, 0, 0, 0};
    vector<double> tempP = {0, 0, 0, 0};
    for (int i = 0; i < schema2; i++) {
        N[i][0] = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
        N[i][1] = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
        N[i][2] = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
        N[i][3] = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
        dNdKsi[i][0] = -0.25 * (1 - eta[i]);
        dNdKsi[i][1] = 0.25 * (1 - eta[i]);
        dNdKsi[i][2] = 0.25 * (1 + eta[i]);
        dNdKsi[i][3] = -0.25 * (1 + eta[i]);
        dNdEta[i][0] = -0.25 * (1 - ksi[i]);
        dNdEta[i][1] = -0.25 * (1 + ksi[i]);
        dNdEta[i][2] = 0.25 * (1 + ksi[i]);
        dNdEta[i][3] = 0.25 * (1 - ksi[i]);
    }
    for (int i = 0; i < schema2; i++) {
        for (int j = 0; j < 4; j++) {
            jacobian[i][0] += dNdKsi[i][j] * b.nodes[j].x;//dxdksi
            jacobian[i][2] += dNdEta[i][j] * b.nodes[j].x;//dxdeta
            jacobian[i][1] += dNdKsi[i][j] * b.nodes[j].y;//dydksi
            jacobian[i][3] += dNdEta[i][j] * b.nodes[j].y;//dydeta
        }
    }
    for (int i = 0; i < schema2; i++) {
        det[i] = (jacobian[i][0] * jacobian[i][3]) - (jacobian[i][1] * jacobian[i][2]);
        reversedJacobian[i][0] = 1 / det[i] * jacobian[i][3];
        reversedJacobian[i][1] = (-1) / det[i] * jacobian[i][1];
        reversedJacobian[i][2] = (-1) / det[i] * jacobian[i][2];
        reversedJacobian[i][3] = 1 / det[i] * jacobian[i][0];
    }
    /*cout << "Element " << b.elemID << " \nDet:" << det[0] << "\tJakobian:\n";
    for (int i = 0; i < 4; i++) {
        cout << jacobian[0][i] << "\t";
    }
    cout << endl << endl;*/

    for (int i = 0; i < schema2; i++) {
        for (int j = 0; j < 4; j++) {
            dNdX[i][j] = (dNdKsi[i][j] * reversedJacobian[i][0] + dNdEta[i][j] * reversedJacobian[i][1]);//git
            dNdY[i][j] = (dNdKsi[i][j] * reversedJacobian[i][2] + dNdEta[i][j] * reversedJacobian[i][3]);//git
        }
    }
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < schema2; b++) {
            dNdXT[a][b] = dNdX[b][a];
            dNdYT[a][b] = dNdY[b][a];
            NT[a][b] = N[b][a];
        }
    }
    for (int point = 1; point <= schema2; point++) {
        /*for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                multipliedX[i][j] = dNdX[point - 1][j] * dNdXT[i][point - 1];
                multipliedY[i][j] = dNdY[point - 1][j] * dNdYT[i][point - 1];
            }
        }*/
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                tempH[i][j] +=
                        grid.heatConductionIndex * det[0] *
                        ((dNdX[point - 1][i] * dNdX[point - 1][j]) + dNdY[point - 1][i] * dNdY[point - 1][j]) *
                        weights[point - 1];//(multipliedX[i][j] + multipliedY[i][j])
            }
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                tempC[i][j] +=
                        det[0] * grid.denisity * grid.c * (N[point - 1][j] * NT[i][point - 1]) *
                        weights[point - 1]; //* multipliedNNT[i][j];
            }
        }
    }
    b.H = sumVectors(b.H, tempH);
    b.C = sumVectors(b.C, tempC);
    //displayArray(tempH);
    //displayArray(tempC);
    vector<vector<double>> Nlocal(schema2, vector<double>(4, 0));
    vector<double> weights2;

    if (grid.schema == 3) {
        weights2 = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    }
    else if (grid.schema == 4) {
        weights2 = {w1,w2,w2,w1};
    }
    /*else if(grid.schema==2){
        weights2={1,1};
    }*/

    //double temp = 1 / (sqrt(3));
    if (b.nodes[0].bc == true && b.nodes[1].bc == true) {
        //pc3 [0][0v1],[1][0v1],[2][0v1]
        for (int i = 0; i < grid.schema; i++) {
            Nlocal[i][0] = 0.5 * (1 - ksi[i]);
            Nlocal[i][1] = 0.5 * (1 + ksi[i]);
        }
        //displayArray(Nlocal,9,4);

        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    HBc[j][k] += ((Nlocal[i][j] * Nlocal[i][k]) * weights2[i]);
                }
            }
        }
        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                tempP[j] +=
                        (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha *
                        weights2[i] * Nlocal[i][j];
            }
        }
        Nlocal.clear();
        Nlocal=vector<vector<double>>(schema2, vector<double>(4, 0));
    }
    if (b.nodes[1].bc == true && b.nodes[2].bc == true) {
        //pc3 [2][1v2],[5][1v2],[8][1v2]

        for (int i = 0; i < grid.schema; i++) {
            Nlocal[grid.schema - 1 + i * grid.schema][1] = 0.5 * (1 - eta[grid.schema - 1 + i * grid.schema]);
            Nlocal[grid.schema - 1 + i * grid.schema][2] = 0.5 * (1 + eta[grid.schema - 1 + i * grid.schema]);
        }
        //displayArray(Nlocal,9,4);
        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    HBc[j][k] += ((Nlocal[grid.schema - 1 + i * grid.schema][j] *
                                   Nlocal[grid.schema - 1 + i * grid.schema][k]) * weights2[i]);
                }
            }
        }
        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                tempP[j] +=
                        (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha *
                        weights2[i] * Nlocal[grid.schema - 1 + i * grid.schema][j];
            }
        }
        Nlocal.clear();
        Nlocal=vector<vector<double>>(schema2, vector<double>(4, 0));
    }
    if (b.nodes[2].bc == true && b.nodes[3].bc == true) {
        //pc3 [6][2v3],[7][2v3],[8][2v3]
        //12,13,14,15
        for (int i = 0; i < grid.schema; i++) {
            Nlocal[schema2-grid.schema + i][2] = 0.5 * (1 + ksi[schema2-grid.schema + i]);
            Nlocal[schema2-grid.schema + i][3] = 0.5 * (1 - ksi[schema2-grid.schema + i]);
        }
        //displayArray(Nlocal,9,4);

        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    HBc[j][k] += ((Nlocal[schema2-grid.schema + i][j] *
                                   Nlocal[schema2-grid.schema + i][k]) * weights2[i]);
                }
            }
        }
        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                tempP[j] +=
                        (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha *
                        weights2[i] * Nlocal[schema2-grid.schema + i][j];

            }
        }
        Nlocal.clear();
        Nlocal=vector<vector<double>>(schema2, vector<double>(4, 0));
    }
    if (b.nodes[3].bc == true && b.nodes[0].bc == true) {
        //pc3 [0][0v3],[3][0v3],[6][0v3]
        for (int i = 0; i < grid.schema; i++) {
            Nlocal[grid.schema * i][0] = 0.5 * (1 - eta[grid.schema * i]);
            Nlocal[grid.schema * i][3] = 0.5 * (1 + eta[grid.schema * i]);
        }
        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    HBc[j][k] += ((Nlocal[grid.schema * i][j] * Nlocal[grid.schema * i][k]) * weights2[i]);
                }
            }
        }
        for (int i = 0; i < grid.schema; i++) {
            for (int j = 0; j < 4; j++) {
                tempP[j] +=
                        (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha *
                        weights2[i] * Nlocal[grid.schema * i][j];
            }
        }
        Nlocal.clear();
        Nlocal=vector<vector<double>>(schema2, vector<double>(4, 0));
    }
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 4; k++) {
            HBc[i][k] *= (grid.alpha * (0.5 * grid.H / (grid.nH - 1)));
            b.H[i][k] += HBc[i][k];
        }
    }
    //cout << "Element " << b.elemID << endl;

    /*cout << "dNdKsi:\n";
    displayArray(dNdKsi, schema2,4);
    cout << "dNdEta:\n";
    displayArray(dNdEta, schema2,4);
    cout << "dNdX:\n";
    displayArray(dNdX, schema2,4);
    cout << "dNdXT:\n";
    displayArray(dNdXT, 4,schema2);
    cout << "dNdY:\n";
    displayArray(dNdY, schema2,4);
    cout << "dNdYT:\n";
    displayArray(dNdYT, 4,schema2);
    cout << "Macierz H+Hbc:\n";
    displayArray(b.H, 4);
    cout << "Macierz C:\n";
    displayArray(b.C, 4);
    cout << "HBC\n";
    displayArray(HBc);
    cout << "P lokalne:\n";
    displayVector(tempP);*/
    LocalMatrixElemData localMatrixElemData;
    for (int i = 0; i < 4; i++) {
        localMatrixElemData.P[i] = tempP[i];
        for (int j = 0; j < 4; j++) {
            localMatrixElemData.H[i][j] = b.H[i][j];
            localMatrixElemData.C[i][j] = b.C[i][j];
        }
    }
    return localMatrixElemData;
}

