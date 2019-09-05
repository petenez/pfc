#!/bin/bash

# This script demonstrates setting up and relaxing a graphene--hexagonal boron nitride heterostructure using PFC and visualizing it

# compile C codes
mpicc pfc-relax.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o pfc-relax
mpicc smooth.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o smooth

# instructions can be printed out by running without input arguments:
# mpirun -np 1 code
# java -jar code.jar

# construct initial state
java -jar pfc-init.jar hetero.i0
java -jar pfc-init.jar hetero.i1
java -jar pfc-init.jar hetero.i2

# relax the system (see the input file hetero.in for details)
mpirun -np 4 pfc-relax hetero.in

for (( i = 0 ; i <= 100000 ; i += 100 )) ; do

	# if data file doesn't exist then skip
	if [[ ! -e hetero-$i.n0 ]] ; then
		continue
	fi
	
	# plot the density fields
	java -jar plotter.jar hetero-$i.n0 hetero-$i-n0.png
	java -jar plotter.jar hetero-$i.n1 hetero-$i-n1.png
	java -jar plotter.jar hetero-$i.n2 hetero-$i-n2.png
	
	# smooth the density fields
	mpirun -np 4 smooth hetero-$i.n0 hetero-$i.eta0 0.2
	mpirun -np 4 smooth hetero-$i.n1 hetero-$i.eta1 0.2
	mpirun -np 4 smooth hetero-$i.n2 hetero-$i.eta2 0.2
	
	# plot the smoothed density fields
	java -jar plotter.jar hetero-$i.eta0 hetero-$i-eta0.png
	java -jar plotter.jar hetero-$i.eta1 hetero-$i-eta1.png
	java -jar plotter.jar hetero-$i.eta2 hetero-$i-eta2.png
	
	# visualize the heterostructure
	java -jar heteroplotter.jar hetero-$i.n0 hetero-$i.n1 hetero-$i.n2 hetero-$i.eta0 hetero-$i.eta1 hetero-$i.eta2 hetero-$i.png
	
	# form sum of density fields
	java -jar manipulator.jar add hetero-$i.n0 hetero-$i.n1 hetero-$i.sum_
	java -jar manipulator.jar add hetero-$i.sum_ hetero-$i.n2 hetero-$i.sum
	rm hetero-$i.sum_
	
	# plot sum
	java -jar plotter.jar hetero-$i.sum hetero-$i-sum.png
done
