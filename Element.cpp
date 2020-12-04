//
// Created by Lasiuk on 04.12.2020.
//

#include "Element.h"
Element::Element() {
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

Element::Element(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
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
Element::Element(Node a, Node b, Node c, Node d) {
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
void Element::displayElement() {
    cout << "\nElement: " << this->elemID << "\tPunkty calkowania:\n";
    for (int i = 0; i < 4; i++) {
        cout.precision(1);
        cout << "Node id:" << this->nodes[i].id + 1 << "  x: " << this->nodes[i].x << "\t\ty: " << this->nodes[i].y<< endl;
        for(int j=0;j<4;j++)
        {
            cout<<"BC flag:"<<this->nodes[j].bc<<"\t";
        }
        cout<<endl;
    }

}
int Element::staticElemID=0;