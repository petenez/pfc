#!/bin/bash

# This script demonstrates setting up and relaxing a PFC system, visualizing it and converting it into atomic coordinates.

# compile C codes
mpicc pfc-relax.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o pfc-relax
mpicc smooth.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o smooth
mpicc fields.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o fields

# instructions can be printed out by running without input arguments:
# mpirun -np 1 code
# java -jar code.jar

# alternatives for constructing the initial state
## initialize with random noise (see the input file honey.i for details)
#java -jar pfc-init.jar honey.i

## initialize with a density field from an image
java -jar img2pfc.jar cameraman.png honey.n
java -jar manipulator.jar r honey.n honey.n 180

# relax the system (see the input file honey.in for details)
mpirun -np 4 pfc-relax honey.in

for (( i = 0 ; i <= 10000 ; i += 100 )) ; do
	
	# plot the density field
	if [[ -e honey-$i.n ]] ; then
		java -jar plotter.jar honey-$i.n honey-$i-n.png
	fi
	
	# smooth the density field
	if [[ -e honey-$i.n ]] ; then
		mpirun -np 4 smooth honey-$i.n honey-$i.eta 0.2
	fi
	
	# plot the smoothed density field
	if [[ -e honey-$i.eta ]] ; then
		java -jar plotter.jar honey-$i.eta honey-$i-eta.png
	fi
	
	# plot the free energy density
	if [[ -e honey-$i.f ]] ; then
		java -jar plotter.jar honey-$i.f honey-$i-f.png
	fi
	
	# smooth the free energy density
	if [[ -e honey-$i.f ]] ; then
		mpirun -np 4 smooth honey-$i.f honey-$i.g 0.2
	fi
	
	# plot the smoothed free energy density
	if [[ -e honey-$i.f ]] ; then
		java -jar plotter.jar honey-$i.g honey-$i-g.png
	fi
	
	
	# get the orientation and deformation fields
	if [[ -e honey-$i.n ]] ; then
		mpirun -np 4 fields honey-$i.n honey-$i.phi honey-$i.chi 6 0.2 0.15 2.0
	fi
	
	# plot the orientation field
	if [[ -e honey-$i.phi ]] ; then
		java -jar plotter.jar honey-$i.phi honey-$i-phi.png
	fi
	
	# plot the deformation field
	if [[ -e honey-$i.chi ]] ; then
		java -jar plotter.jar honey-$i.chi honey-$i-chi.png
	fi
	
	# plot the oriented density field
	if [[ -e honey-$i.n && -e honey-$i.phi ]] ; then
		java -jar orienter.jar honey-$i.n honey-$i.phi honey-$i-nphi.png
	fi
	
done
	
if [[ -e honey-10000.n ]] ; then

	# convert the density field into atomic coordinates
	java -jar coordinator.jar honey-10000.n honey-10000.xyz 7.26 2.46

	# create a style file to color the atoms based on their bond angles
	java -jar colorizer.jar honey-10000.xyz honey-10000.sty 0.5 0.15

	# create a POV-Ray script based on the coordinates
	java -jar pover.jar honey.par

	# render using POV-Ray (latter with transparency (if specified by the user))
	povray +ihoney-10000.pov +ohoney-10000-pov.png +w1000 +h1000 -d
	#povray +ihoney-10000.pov +ohoney-10000-pov.png +w1000 +h1000 -d +ua

fi
