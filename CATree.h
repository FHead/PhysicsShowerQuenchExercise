#include <vector>
#include <map>

#include "TauHelperFunctions3.h"

#define TYPE_QUARK 0
#define TYPE_GLUON 1

class Node;
struct NodePair;
void BuildCATree(std::vector<Node *> &Nodes, double p);
void BuildCATree2(std::vector<Node *> &Nodes);
NodePair FindClosestPair(std::vector<Node *> &Nodes, std::vector<std::pair<double, int>> &NodeEta);
Node *FindSDNode(Node *HeadNode, double ZCut = 0.1, double Beta = 0, double R0 = 0.4);
std::vector<std::pair<double, double>> CountSD(Node *HeadNode, double ZCut = 0.1, double Beta = 0, double R0 = 0.4, double DRCut = 0.1);
int NodeDistance(Node *Child, Node *Root);
double SDCSum(std::vector<std::pair<double, double>> &Z, double Kappa);
std::pair<double, double> WinnerTakesAllAxis(Node *N);
void AssignTime(Node *N);
void AssignQG(Node *N);

class Node
{
public:
   FourVector P;
   FourVector V;
   int QG;
   Node *Child1;
   Node *Child2;
   Node *Parent;
   int N;   // multiplicity: how many nodes (including self) under this node
public:
   Node();
   Node(FourVector &p);
   Node(Node *n1, Node *n2);
   ~Node();
};

struct NodePair
{
   int ID1;
   int ID2;
   Node *N1;
   Node *N2;
   double DR2;
   double DEta;
   double DPhi;
};



