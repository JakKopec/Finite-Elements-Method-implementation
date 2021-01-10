#include <cmath>
#include "Elem2Solve.h"
LocalMatrixElem2 elem2solve(Element b,FEMGrid grid) {
    cout.precision(3);
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

    vector<vector<double>> HBc = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };
    vector<double> NBc1 = {0, 0, 0, 0};
    vector<double> NBc2 = {0, 0, 0, 0};

    vector<double> tempP = {0, 0, 0, 0};

    for (int i = 0; i < 4; i++) {
        jacobian[0] += dNdKsi[0][i] * b.nodes[i].x;
        jacobian[1] += dNdEta[1][i] * b.nodes[i].x;
        jacobian[2] += dNdKsi[2][i] * b.nodes[i].y;
        jacobian[3] += dNdEta[3][i] * b.nodes[i].y;
    }
    det = (jacobian[0] * jacobian[3]) - (jacobian[1] * jacobian[2]);
    reversedJacobian[0] = jacobian[3] / det;
    reversedJacobian[1] = -jacobian[1] / det;
    reversedJacobian[2] = -jacobian[2] / det;
    reversedJacobian[3] = jacobian[0] / det;
    /*cout << "Element " << b.elemID << " \nDet:" << det << "\tJakobian:\n";
    for (int j = 0; j < 4; j++) {
        cout << jacobian[j] << "\t";
    }
    cout << endl;*/

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
                tempH[i][j] =
                        grid.heatConductionIndex * det * (multipliedX[i][j] + multipliedY[i][j]) * weight[point - 1];
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                tempC[i][j] =
                        det * grid.denisity * grid.c * (N[point - 1][j] * NT[point - 1][i]); //* multipliedNNT[i][j];

            }
        }
        b.H = sumVectors(b.H, tempH);
        b.C = sumVectors(b.C, tempC);
        tempH = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
        tempC = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
    }

    vector<vector<double>> Nlocal = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    };


    double temp = 1 / (sqrt(3));
    if (b.nodes[0].bc == true && b.nodes[1].bc == true) {
        Nlocal[0][0] = shapeFun(-temp, -1, -1, -1);;
        Nlocal[0][1] = shapeFun(-temp, -1, 1, -1);
        Nlocal[1][0] = shapeFun(temp, -1, -1, -1);
        Nlocal[1][1] = shapeFun(temp, -1, 1, -1);
        for (int a = 0; a < 2; a++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBc[i][j] += (Nlocal[a][i] * Nlocal[a][j]);
                }
            }
        }
        vector<double> vec;
        vec.push_back(Nlocal[0][0]);
        vec.push_back(Nlocal[0][1]);
        vec.push_back(Nlocal[1][0]);
        vec.push_back(Nlocal[1][1]);
        tempP[0] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[0] + vec[1]);
        tempP[1] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[2] + vec[3]);
        vec.clear();
        Nlocal = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };


    }

    if (b.nodes[1].bc == true && b.nodes[2].bc == true) {
        Nlocal[1][1] = shapeFun(1, -temp, 1, -1);
        Nlocal[1][2] = shapeFun(1, -temp, 1, 1);
        Nlocal[2][1] = shapeFun(1, temp, 1, -1);
        Nlocal[2][2] = shapeFun(1, temp, 1, 1);
        for (int a = 0; a < 2; a++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBc[i][j] += (Nlocal[1 + a][i] * Nlocal[1 + a][j]);
                }
            }
        }
        vector<double> vec;
        vec.push_back(Nlocal[1][1]);
        vec.push_back(Nlocal[1][2]);
        vec.push_back(Nlocal[2][1]);
        vec.push_back(Nlocal[2][2]);
        tempP[1] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[0] + vec[1]);
        tempP[2] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[2] + vec[3]);
        vec.clear();
        Nlocal = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
    }
    if (b.nodes[2].bc == true && b.nodes[3].bc == true) {
        Nlocal[2][2] = shapeFun(temp, 1, 1, 1);
        Nlocal[2][3] = shapeFun(temp, 1, -1, 1);
        Nlocal[3][2] = shapeFun(-temp, 1, 1, 1);
        Nlocal[3][3] = shapeFun(-temp, 1, -1, 1);
        for (int a = 0; a < 2; a++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBc[i][j] += (Nlocal[2 + a][i] * Nlocal[2 + a][j]);
                }
            }
        }
        vector<double> vec;
        vec.push_back(Nlocal[2][2]);
        vec.push_back(Nlocal[2][3]);
        vec.push_back(Nlocal[3][2]);
        vec.push_back(Nlocal[3][3]);
        tempP[2] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[0] + vec[1]);
        tempP[3] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[2] + vec[3]);
        vec.clear();
        Nlocal = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
    }
    if (b.nodes[3].bc == true && b.nodes[0].bc == true) {
        Nlocal[0][0] = shapeFun(-1, temp, -1, 1);//0.5 * (1 - eta[0]);
        Nlocal[0][3] = shapeFun(-1, temp, -1, -1);//0.5 * (1 + eta[3]);
        Nlocal[3][0] = shapeFun(-1, -temp, -1, 1);
        Nlocal[3][3] = shapeFun(-1, -temp, -1, -1);
        for (int a = 0; a < 2; a++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBc[i][j] += (Nlocal[3 * a][i] * Nlocal[3 * a][j]);
                }
            }
        }
        vector<double> vec;
        vec.push_back(Nlocal[0][0]);
        vec.push_back(Nlocal[0][3]);
        vec.push_back(Nlocal[3][0]);
        vec.push_back(Nlocal[3][3]);
        tempP[0] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[0] + vec[1]);
        tempP[3] +=
                 (0.5 * grid.H / (grid.nH - 1)) * grid.ambientTemperature * grid.alpha * (vec[2] + vec[3]);
        vec.clear();
        Nlocal = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            HBc[i][j] *= (grid.alpha * (0.5 * grid.H / (grid.nH - 1)));
            b.H[i][j] += HBc[i][j];
        }
    }

    /*cout << "Element " << b.elemID << endl;

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
    displayArray(multipliedY,4);
    cout << "Macierz H+Hbc:\n";
    displayArray(b.H, 4);
    cout << "Macierz C:\n";
    displayArray(b.C, 4)
    cout << "HBC\n";
    displayArray(HBc);
    cout << "P lokalne:\n";
    displayVector(tempP);*/

    LocalMatrixElem2 localMatrixElem2;
    for (int i = 0; i < 4; i++) {
        localMatrixElem2.P[i] = tempP[i];
        for (int j = 0; j < 4; j++) {
            localMatrixElem2.H[i][j] = b.H[i][j];
            localMatrixElem2.C[i][j] = b.C[i][j];
        }
    }

    return localMatrixElem2;
}
LocalMatrixElem2 elem3solve(Element b,FEMGrid grid) {
    double sq = sqrt(15) / 5;
    double ksi[9] =
            {-sq, 0, sq, -sq, 0, sq, -sq, 0, sq};//do poprawy
    double eta[9] =
            {-sq, -sq, -sq, 0, 0, 0, sq, sq, sq};//do poprawy
    int weight[9] =
            {25 / 81, 40 / 81, 25 / 81, 40 / 81, 64 / 81, 40 / 81, 25 / 81, 40 / 81, 25 / 81};
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
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 4; j++) {
            jacobian[i][0] += dNdKsi[0][j] * b.nodes[j].x;
            jacobian[i][1] += dNdEta[1][j] * b.nodes[j].x;
            jacobian[i][2] += dNdKsi[2][j] * b.nodes[j].y;
            jacobian[i][3] += dNdEta[3][j] * b.nodes[j].y;
        }

        det[i] = (jacobian[i][0] * jacobian[i][3]) - (jacobian[i][2] * jacobian[i][1]);
        reversedJacobian[i][0] = jacobian[i][3] / det[i];
        reversedJacobian[i][1] = (-1)*jacobian[i][1] / det[i];
        reversedJacobian[i][2] = (-1)*jacobian[i][2] / det[i];
        reversedJacobian[i][3] = jacobian[i][0] / det[i];
    }
    cout << "Jacobian\n";
    displayArray(jacobian);
    cout << "Reversed jacobian\n";
    displayArray(reversedJacobian);
    cout << "Det\n";
    for (int i = 0; i < 9; i++) cout << det[i] << endl;
    vector<vector<double>> dNdX = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    vector<vector<double>> dNdY = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    for (int a = 0; a < 9; a++) {
        cout << "dNdKsi:\n";
        displayArray(dNdKsi, 4);
        cout << "dNdEta:\n";
        displayArray(dNdEta, 4);
        dNdX[0][a] = (dNdKsi[0][a] * reversedJacobian[a][0] + dNdEta[0][a] * reversedJacobian[a][1]);
        dNdX[1][a] = (dNdKsi[1][a] * reversedJacobian[a][0] + dNdEta[1][a] * reversedJacobian[a][1]);
        dNdX[2][a] = (dNdKsi[2][a] * reversedJacobian[a][0] + dNdEta[2][a] * reversedJacobian[a][1]);
        dNdX[3][a] = (dNdKsi[3][a] * reversedJacobian[a][0] + dNdEta[3][a] * reversedJacobian[a][1]);

        dNdY[0][a] = (dNdKsi[0][a] * reversedJacobian[a][2] + dNdEta[0][a] * reversedJacobian[a][3]);
        dNdY[1][a] = (dNdKsi[1][a] * reversedJacobian[a][2] + dNdEta[1][a] * reversedJacobian[a][3]);
        dNdY[2][a] = (dNdKsi[2][a] * reversedJacobian[a][2] + dNdEta[2][a] * reversedJacobian[a][3]);
        dNdY[3][a] = (dNdKsi[3][a] * reversedJacobian[a][2] + dNdEta[3][a] * reversedJacobian[a][3]);
    }
    vector<vector<double>> dNdXT = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<vector<double>> dNdYT = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    for (int a = 0; a < 9; a++) {
        for (int b = 0; b < 4; b++) {
            dNdXT[a][b] = dNdX[b][a];
            dNdYT[a][b] = dNdY[b][a];
        }
    }
    vector<vector<double>> N = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<vector<double>> NT = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    for (int i = 0; i < 9; i++) {
        N[i][0] = (0.25 * (1 - ksi[i]) * (1 - eta[i]));
        N[i][1] = (0.25 * (1 + ksi[i]) * (1 - eta[i]));
        N[i][2] = (0.25 * (1 + ksi[i]) * (1 + eta[i]));
        N[i][3] = (0.25 * (1 - ksi[i]) * (1 + eta[i]));
    }
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 9; b++) {
            NT[a][b] = N[b][a];
        }
    }
    vector<vector<double>> multipliedX = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    vector<vector<double>> multipliedY = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 4; k++) {
                multipliedX[i][j] += dNdX[k][i] * dNdXT[j][k];
                multipliedY[i][j] += dNdY[k][i] * dNdYT[j][k];
            }
        }
    }
    //displayArray(multipliedX);
    cout<<endl<<endl;
    //displayArray(multipliedY);
    vector<vector<double>> tempH = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<vector<double>> tempC = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    vector<vector<double>> C = {
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    for (int point = 1; point <= 9; point++) {
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 9; j++) {
                multipliedX[i][j] = dNdX[point - 1][j] * dNdXT[i][point - 1];
                multipliedY[i][j] = dNdY[point - 1][j] * dNdYT[i][point - 1];
            }
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                tempH[i][j] = grid.heatConductionIndex * det[point - 1] * (multipliedX[i][j] + multipliedY[i][j]) *
                              weight[point - 1];
            }
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                tempC[i][j] = det[point - 1] * grid.denisity * grid.c * (N[point - 1][j] * NT[i][point - 1]);
                //cout << tempC[i][j] << "\t";
            }
            //cout << endl;
        }
        cout<<"C\n";
        displayArray(tempC);
        cout<<endl;
        cout<<"H\n";
        displayArray(tempH);
    }



/*
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
//}
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
displayArray(multipliedY,4);
cout << "Macierz H:\n";
displayArray(b.H);
cout << "Macierz C:\n";
displayArray(b.C);
*/

    LocalMatrixElem2 localMatrixElem2;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            localMatrixElem2.H[i][j] = b.H[i][j];
            localMatrixElem2.C[i][j] = b.C[i][j];
        }
    }
    return localMatrixElem2;
}