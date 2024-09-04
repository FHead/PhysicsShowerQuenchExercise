#include <iostream>
using namespace std;

#include "TauHelperFunctions3.h"
#include "CATree.h"
#include "CommandLine.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   vector<FourVector> Particles;

   double P[] = {0, 0.5, 1};

   for(double p : P)
   {
      cout << p << endl;

      vector<Node *> Nodes;
      for(FourVector &P : Particles)
         Nodes.push_back(new Node(P));

      BuildCATree(Nodes, p);

      AssignTime(Nodes[0]);

      delete Nodes[0];
   }

   return 0;
}



