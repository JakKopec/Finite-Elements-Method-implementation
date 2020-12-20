#ifndef UNTITLED_ELEMENT_H
#define UNTITLED_ELEMENT_H
#include "Node.h"
#include<vector>
#include <iostream>

using namespace std;

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
    vector<double> P;
    Element();
    Element(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4);
    Element(Node a, Node b, Node c, Node d);
    void displayElement();
};


#endif //UNTITLED_ELEMENT_H
