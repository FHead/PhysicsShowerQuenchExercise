
# Energy loss integration on a pre-existing shower structure

## Introduction

We start from some preexisting hepmc2 input files, and the code runs
- fastjet clustering with anti-kt R = 0.4
- for each jet we take all the constitents and connect shower history with generalized kt p = 0.0, 0.5, 1.0
- for each setting, integrate some pre-existing dE/dx on the lines
- we report the outputs in the test.out text file


## Practically

Clone this repository and do ``make''.  You will need fastjet installed (if fastjet-config is in the path it's ok)


## Output file format

7 columns: jet PT, jet eta, jet phi, constituents, eloss p = 0, eloss p = 0.5, eloss p = 1
each row is a jet


## Some more information







