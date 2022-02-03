#!/bin/bash

# begin post-processing
ls -d plt* >plotfiles.txt
declare -i cores=4
#Iterate over the plot files 
while IFS= read -r pltfile; do
	echo "${r}=================================================${nc}"
        echo "${g}PROCESSING ${y}$pltfile...${nc}"
        echo "${r}Computing transport coeffs, thermal diffusivity and mixture fraction..${nc}"
        mpirun -np $cores ./plotTransportAlphaZ2d.gnu.MPI.ex infile=$pltfile < /dev/null
        echo "${g}Done computing transport ${nc}"
        echo "${r}Appending transport to plot file..${nc}"
        b="_D"
        d="_sdr"
        c=$pltfile
        c+=$b
        e=$pltfile
        e+=$d
        echo $c
        mpirun -np $cores ./combinePlts2d.gnu.MPI.ex infileR=$c infileL=$pltfile outfile =$pltfile < /dev/null
        rm -r plt*_D
        echo "${g}Done dumping data in plot file${nc}"
        echo "${r}Computing scalar dissipation rate..${nc}"
        mpirun -np $cores ./scalar_dissipation_rate2d.gnu.MPI.ex inp.sdr infile=$pltfile outfile=$e < /dev/null
        echo "${g}Done computing scalar dissipation rate ${nc}"
        echo "${r}Appending scalar dissipation rate to plot file..${nc}"
        mpirun -np $cores ./combinePlts2d.gnu.MPI.ex infileR=$e infileL=$pltfile outfile =$pltfile < /dev/null
        rm -r plt*_sdr
        echo "${g}Done dumping data in plot file${nc}"
        echo "${r}Computing conditional means..${nc}"
        mpirun -np $cores ./conditionalMean2d.gnu.MPI.ex inp.cond_mean infile=$pltfile < /dev/null
        echo "${g}Done Computing conditional means..${nc}"
	echo "${r}Writing data for scatter plots..${nc}"
        ./scatterPlot2d.gnu.MPI.ex inp.scatter Plot_File=$pltfile < /dev/null
	mv scatter.dat $pltfile
	echo "${g}Done writing scatter data...${nc}"
	echo "${r}Time-tagging scatter and conditional plot data"
	cd $pltfile
	e="_scatter.dat"
	f=$pltfile
	f+=$e
	mv scatter.dat $f
        q="_CM.dat"	
	h=$pltfile
	h+=$q
	mv CM*dat $h
	mv CM*key CM.key
	cd ../
        echo -e "${r}Done processing $pltfile \n....${nc}"
        echo "${r}=================================================${nc}"
done < plotfiles.txt
echo "${r}POST PROCESSING COMPLETE..${nc}"
echo "${r}=================================================${nc}"
echo "${g} Please use the MATLAB scripts (.m) files to generate the scatter plots${nc}"
echo "${r}=================================================${nc}"
