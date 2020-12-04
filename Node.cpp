#include <iostream>
#include "Node.h"
using namespace std;
void Node::Multiply(double dw,double dh,int muldw,int muldh)//przypisanie wspolrzednych wezlom
{
    this->x=(dw*muldw);
    this->y=(dh*muldh);
}
Node& Node::operator=(const Node &a)
{
    this->x=a.x;
    this->y=a.y;
    this->id=a.id;
    this->bc=a.bc;
    return *this;
}
void Node::displayNode()
{
    cout.precision(1);
    cout<<"Node ID:"<<id+1<<"\tx="<<this->x<<"\t\ty="<<this->y<<"\t\tt0="<<this->t0<<"\t\tflaga wb:"<<bc<<endl;
}
int Node::returnId(){
    return id;
}
Node::Node(double a,double b)
{
    this->x=a;
    this->y=b;
}
Node::Node(double a,double b,int id)
{
    this->x=a;
    this->y=b;
    this->id=id;
}