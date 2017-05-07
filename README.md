# SumProductGraphMod - Quick Start

Hello! This Library implements Sum-Product Graphical Models (SPGM), including inference and structure learning.

This document contains instructions on how to run the code and reproduce the experiments in the paper. 
Detailed examples can be found in the Test units alongside each C++ class. 

## 1) running / compiling the code

If you don't need to compile the code, you can just use the executable located in dist/Release/GNU-Linux-x86. 

You can compile the C++ project with either Netbeans 8.2 (the folder is a Netbeans project) or with Make (run make in the main folder).
The only library you need is the BOOST C++ library ( http://www.boost.org/ ). 


## 2) reproducing the experiments in the paper

The program executing LearnSPGM is called spgm_0.2 and can be found in dist/Release/GNU-Linux-x86. 

To run the experiments in the paper you have first to download the data from http://spn.cs.washington.edu/learnspn/ 

then go to the directory where spgm_0.2 is located and run the following command:

./spgm_0.2  --dataFolder YOUR_DATA_FOLDER_PATH --reproduce  

the last argument reproduces learning with the hyperparameters used in the paper. 

## 3) running your own experiments

spgm_0.2 accepts a number of parameters (type --help to see a dirty list of them). The most important are:
--nIns (number of edge insertions)
--K (number of mixture components in mixtures of SPGMs)
--dir (dirichelet prior)
--wp (weight prior)
--name (dataset name, chosen from "nltcs", "msnbc", "kdd", "plants", "baudio", "jester", "bnetflix", "accidents", "tretail", "pumsb_star", "dna", "kosarek", "msweb", "book", "tmovie", "cwebkb", "cr52", "c20ng", "bbc", "ad" )

Enjoy!
