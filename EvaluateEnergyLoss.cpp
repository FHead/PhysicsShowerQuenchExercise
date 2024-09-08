#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

#include "TauHelperFunctions3.h"
#include "CATree.h"
#include "CommandLine.h"

int main(int argc, char *argv[]);
void RunAnalysis(vector<PseudoJet> &V, ofstream &out);
double dEdx(int Type, double T, double x);
double ParticleELoss(Node *N);
double GetT(double Time);
double GetTotalELoss(Node *N);

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName = CL.Get("Input");
   string OutputFileName = CL.Get("Output", "test.out");

   ofstream out(OutputFileName);

   vector<PseudoJet> V;

   ifstream in(InputFileName);
   while(in)
   {
      char ch[10240];
      ch[0] = '\0';
      in.getline(ch, 10239, '\n');
      if(ch[0] == '\0')
         continue;

      stringstream str(ch);

      vector<string> list;
      while(str)
      {
         string Temp = "";
         str >> Temp;
         if(Temp != "")
            list.push_back(Temp);
      }

      if(list.size() == 0)
         continue;
      if(list[0][0] == '#')
         continue;
      if(list[0].size() > 1)
         continue;
      if(list[0] == "E")
      {
         if(V.size() > 0)
            RunAnalysis(V, out);
         V.clear();
         // MZHadron.EventWeight = stof(list[list.size()-1]);
         continue;
      }
      if(list[0] == "P")
      {
         int ID     = stoi(list[2]);
         double px  = stof(list[3]);
         double py  = stof(list[4]);
         double pz  = stof(list[5]);
         double e   = stof(list[6]);
         double m   = stof(list[7]);
         int Status = stoi(list[8]);

         PseudoJet P(px, py, pz, e);
         P.set_user_index(ID);

         if(Status == 1)
            V.push_back(P);
      }
   }
   if(V.size() > 0)
      RunAnalysis(V, out);

   out.close();

   return 0;
}

void RunAnalysis(vector<PseudoJet> &V, ofstream &out)
{
   JetDefinition Definition(antikt_algorithm, 0.4);
   ClusterSequence Sequence(V, Definition);
   vector<PseudoJet> FastJets = sorted_by_pt(Sequence.inclusive_jets(5));

   for(int iJ = 0; iJ < (int)FastJets.size(); iJ++)
   {
      vector<FourVector> Particles;
      vector<bool> IsQuark;
      PseudoJet Jet = FastJets[iJ];
      vector<PseudoJet> Constituents = FastJets[iJ].constituents();
      for(int i = 0; i < Constituents.size(); i++)
      {
         PseudoJet &J = Constituents[i];
         FourVector P(J.e(), J.px(), J.py(), J.pz());
         Particles.push_back(P);
         IsQuark.push_back(abs(J.user_index()) <= 6);
      }

      if(Particles.size() == 0)
         continue;

      // cout << Jet.perp() << endl;

      out << Jet.perp() << " " << Jet.eta() << " " << Jet.phi() << " " << Constituents.size();

      double PValues[] = {0.0, 0.5, 1.0};

      for(double p : PValues)
      {
         // cout << p << endl;

         vector<Node *> Nodes;
         for(int i = 0; i < (int)Particles.size(); i++)
         {
            Node *N = new Node(Particles[i]);
            N->QG = IsQuark[i] ? TYPE_QUARK : TYPE_GLUON;
            Nodes.push_back(N);
         }

         BuildCATree(Nodes, p);
         AssignTime(Nodes[0]);
         AssignQG(Nodes[0]);

         // Now we extract information we want from the tree

         out << " " << GetTotalELoss(Nodes[0]);

         delete Nodes[0];
      }

      out << endl;
   }
}

double dEdx(int Type, double T, double x)
{
   double KRad = 4;
   double CRCF = (Type == TYPE_GLUON) ? 2.25 : 1;

   return KRad * CRCF * T * T * T * x * 0.197 * 0.197;
   
   // double KColl = 2.5;
   
   // return KColl * CRCF * T * T * 0.197 * 0.5;
   
   // if(T > 0)
   //    return 1;
   // return 0;
}

double ParticleELoss(Node *N)
{
   if(N == nullptr)
      return 0;

   double ELoss = 0;

   double E = N->P[0];
   double Q2 = N->P.GetMass2();
   double t = (Q2 > 0) ? (2 * E / Q2 * 0.197) : 50;   // set to some large number (=50) if virtuality is 0
   // t = t * 0.197;

   double XInitial = N->V[0];
   double XFinal = XInitial + t;

   // brute force for now
   double NStep = 2000;
   for(int i = 0; i < NStep; i++)
   {
      double DX = (XFinal - XInitial) / NStep;
      double X = XInitial + DX * (i + 0.5);

      double DE = dEdx(N->QG, GetT(X), X - XInitial) * DX;

      ELoss = ELoss + DE;
   }

   return ELoss;
}

double GetT(double Time)
{
   double FinalTime = 15;   // after this QGP cease to exist, in fm

   if(Time > FinalTime)
      return 0;

   double T0 = 0.50;   // initial temperature
   double T1 = 0.15;   // final temperature before jumping down to 0

   return T0 + (T1 - T0) * Time / FinalTime;
}

double GetTotalELoss(Node *N)
{
   if(N == nullptr)
      return 0;

   return GetTotalELoss(N->Child1) + GetTotalELoss(N->Child2) + ParticleELoss(N);
   // return GetTotalELoss(N->Child1) + ParticleELoss(N);
}



