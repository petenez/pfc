// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

// compile: mpicc smooth.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o smooth

// print instructions: mpirun -np 1 smooth

#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define pi 3.14159265358979323846264338327

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	fftw_mpi_init();
	
	// print instructions
	if(argc == 1) {
		int id;
		MPI_Comm_rank(MPI_COMM_WORLD, &id);
		if(id == 0) {
			printf("This is a tool for smoothing PFC data. The program expects the names of the input and output data files, and the spread of the smoothing kernel in reciprocal units.\n\n");
			printf("An example of valid syntax:\n\n");
			printf("mpirun -np 4 smooth [input data file] [output data file] [smoothing kernel spread]\n");
		 }
		 
		MPI_Finalize();
		return 0;
	}
	
	// input file exists?
	char name[BUFSIZ];
	sprintf(name, "%s", argv[1]);
	FILE *file = fopen(name, "r");
	if(file == NULL) {
		printf("Input file not found!\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	// initialize arrays
	int W, H;
	double dx, dy;
	ptrdiff_t lW, lH, li0, lj0;
	if(fscanf(file, " %d %d %lf %lf ", &W, &H, &dx, &dy));
	int lWHh = fftw_mpi_local_size_2d_transposed(H, W/2 + 1, MPI_COMM_WORLD, &lH, &lj0, &lW, &li0);
	double *f = (double*)fftw_malloc(2*lWHh*sizeof(double));	// data
	double *g = (double*)fftw_malloc(2*lWHh*sizeof(double));	// Gaussian
	int Wp = 2*(W/2 + 1);
	
	// read data
	int i, j, gi, k;
	double dummy;
	for(j = 0; j < H; j++) {
		if(j >= lj0 && j < lj0 + lH) {
			k = Wp*(j - lj0);
			for(i = 0; i < W; i++) {
				if(fscanf(file, " %lf ", &f[k + i]));
			}
		}
		else {
			for(i = 0; i < W; i++) {
				if(fscanf(file, " %lf ", &dummy));
			}
		}
	}
	fclose(file);
	
	// initialize plans
	fftw_plan f_F = fftw_mpi_plan_dft_r2c_2d(H, W, f, (fftw_complex*)f, MPI_COMM_WORLD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_OUT);
	fftw_plan F_f = fftw_mpi_plan_dft_c2r_2d(H, W, (fftw_complex*)f, f, MPI_COMM_WORLD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);
	fftw_plan G_g = fftw_mpi_plan_dft_c2r_2d(H, W, (fftw_complex*)g, g, MPI_COMM_WORLD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);
	
	// convolve in Fourier space
	double WH = W*H;
	double dkx = 2.0*pi/(dx*W);
	double dky = 2.0*pi/(dy*H);
	double sigma = atof(argv[3]);
	double c = 0.5/sigma/sigma;
	double kx2, ky, k2;
	fftw_complex *F = (fftw_complex*)&f[0];
	fftw_complex *G = (fftw_complex*)&g[0];
	fftw_execute(f_F);
	for(i = 0; i < lW; i++) {
		gi = li0 + i;
		kx2 = gi*dkx;
		kx2 *= kx2;
		k = H*i;
		for(j = 0; j < H; j++) {
			if(j < H/2) ky = j*dky;
			else ky = (j - H)*dky;
			k2 = kx2 + ky*ky;
			G[k + j] = exp(- c*k2)/WH;
			F[k + j] *= G[k + j];
		}
	}
	
	// return to direct space
	fftw_execute(F_f);
	fftw_execute(G_g);
	
	// normalize
	double S = 0.0;
	for(j = 0; j < lH; j++) {
		k = Wp*j;
		for(i = 0; i < W; i++) {
			S += g[k + i];
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &S, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(j = 0; j < lH; j++) {
		k = Wp*j;
		for(i = 0; i < W; i++) {
			f[k + i] /= S;
		}
	}
	
	// write output
	sprintf(name, "%s", argv[2]);
	int id, p, P;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	for(p = 0; p < P; p++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if(id == p) {
			if(id == 0) {
				file = fopen(name, "w");
				fprintf(file, "%d %d %lf %lf\n", W, H, dx, dy);
			}
			else file = fopen(name, "a");
			for(j = 0; j < lH; j++) {
				k = Wp*j;
				for(i = 0; i < W; i++) {
					fprintf(file, "%.12e\n", f[k + i]);
				}
			}
			fclose(file);
		}
	}
	
	// clean up
	fftw_free(f);
	fftw_free(g);
	fftw_destroy_plan(f_F);
	fftw_destroy_plan(F_f);
	fftw_destroy_plan(G_g);
	MPI_Finalize();
	return 0;
}
