# an input file for pfc-init
# for instructions, run: java -jar pfc-init.jar

# grid size and discretization
D

#	Nx		Ny		dx		dy
	512		512		0.7		0.7

# Voronoi grains
V	1

# 1st grain
#	x		y		type
	0.0		0.0		0
	
#	sym		a		nave	amp		namp	
	3		7.2552	0.06	0.0		1.0
	
#	tx1		ty1		theta	xr		yr		tx2		ty2
	0.0		0.0		0.0		0.0		0.0		0.0		0.0

# transform
T

#	nave	famp
	0.06	1.0

# output
O

#	filename	precision
	honey.n		6
