#ifndef UNTITLED_NODE_H
#define UNTITLED_NODE_H

class Node
{
public:
    double x;
    double y;
    int id;
    double t0;
    bool bc=false;
    void Multiply(double dw,double dh,int muldw,int muldh);//przypisanie wspolrzednych wezlom
    Node& operator=(const Node &a);
    void displayNode();
    int returnId();
    Node(){};
    Node(double a,double b);
    Node(double a,double b,int id);
};


#endif //UNTITLED_NODE_H
