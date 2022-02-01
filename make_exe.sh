#!/bin/bash

r=`tput setaf 1`
g=`tput setaf 2`
y=`tput setaf 4`
nc=`tput sgr0` 

echo "${r}=================================================${nc}"
echo "${r}Enter number of cores available for compiling${nc}"
read cores
echo "${r}=================================================${nc}"

echo "${r}CLEANING...{nc}"
make clean

echo "${r}=================================================${nc}"
echo "${r}DONE CLEANING...${nc}"
echo "${r}=================================================${nc}"

echo "${r}COMPILING scalar_dissipation_rate.cpp...${nc}"
echo "${r}=================================================${nc}"
make -j$cores EBASE=scalar_dissipation_rate USE_MPI=TRUE
echo "${r}=================================================${nc}"
echo "${r}DONE COMPILING scalar_dissipation_rate.cpp...${nc}"

echo "${r}COMPILING combinePlts.cpp...${nc}"
echo "${r}=================================================${nc}"
make -j$cores EBASE=combinePlts USE_MPI=TRUE
echo "${r}=================================================${nc}"
echo "${r}DONE COMPILING combinePlts.cpp...${nc}"

echo "${r}COMPILING conditionalMean.cpp...${nc}"
echo "${r}=================================================${nc}"
make -j$cores EBASE=conditionalMean USE_MPI=TRUE
echo "${r}=================================================${nc}"
echo "${r}DONE COMPILING conditionalMean.cpp...${nc}"

echo "${r}COMPILING scatterPlot.cpp...${nc}"
echo "${r}=================================================${nc}"
make -j$cores EBASE=scatterPlot USE_MPI=TRUE
echo "${r}=================================================${nc}"
echo "${r}DONE COMPILING scatterPlot.cpp...${nc}"

echo "${r}=================================================${nc}"
echo "${r}COMPILING plotTransportAlphaZ.cpp...${nc}"
cd ModelSpecificAnalysis
make clean
make -j$cores EBASE=plotTransportAlphaZ USE_MPI=TRUE
scp -r *ex ../
cd ../
echo "${r}DONE COMPILING plotTransportAlphaZ.cpp...${nc}"
echo "${r}=================================================${nc}"
echo "${r}COLLECTING EXECUTABLES.....${nc}"
echo "${r}=================================================${nc}"
#make a folder with all the executables and input files
mkdir executables
scp -r *ex executables
echo "${r}=================================================${nc}"
echo "${r}COLLECTING EXECUTABLES.....${nc}"
echo "${r}=================================================${nc}"

#create copies of the input files inside the folder
scp -r inp.sdr inp.scatter inp.cond_mean executables

#move the plotfiles to the folder containing exectuables
mv plt* postProcess.sh executables
cd executables

echo "${r}=================================================${nc}"
echo "${r}READY TO POST-PROCESS.....${nc}"
echo "${r}=================================================${nc}"
