
# Energy loss integration on a pre-existing shower structure

## Introduction

We start from some preexisting hepmc2 input files, and the code runs
- fastjet clustering with anti-kt R = 0.4
- for each jet we take all the constitents and connect shower history with generalized kt p = 0.0, 0.5, 1.0
- for each setting, integrate some pre-existing dE/dx on the lines
- we report the outputs in the test.out text file

If you want to try different new scenarios yourself, we typically only need to change the functions (T profile & dEdx) in the EvaluateEnergyLoss.cpp.  The rest are helper functions.



## Practically

Clone this repository and do `make`.  You will need fastjet installed (if fastjet-config is in the path it is ok)


## Output file format

7 columns: jet PT, jet eta, jet phi, constituents, eloss p = 0, eloss p = 0.5, eloss p = 1
each row is a jet


## Some more information

### CATree.h/CATree.cpp

Here there are a few relevant functions
- BuildCATree(vector<>, double): builds the parton shower tree.  Second parameter is the p in the generalized kt
- AssignTime(...): go through the tree and assign time.  New function for this exercise.  In principle if things look good we don't need to touch this anymore.
- AssignQG(...): go through the tree and assign q/g type.  If the two children are qq or gg, assign g.  If qg or gq, assign q.

### EvaluateEnergyLoss.cpp

The main program that does the energy loss calculation.  Functions...
- main(...): the main function that reads hepmc2 file and calls RunAnalysis(...)
- RunAnalysis(...): cluster jet with anti-kt, then for each clustered jet integrate energy loss
- dEdx(...): calculates dE/dx given some parameters
- ParticleELoss(...): integrate over one internal particle line
- GetT(...): gets the temperature given some time.  Change medium assumptions here
- GetTotalEnergyLoss(...): calls ParticleELoss(...) recursively and adds up everything





