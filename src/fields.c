// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

// compile: mpicc fields.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o fields

// print instructions: mpirun -np 1 fields

#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define pi 3.14159265358979323846264338327

// contains parameters, arrays and FFTW plans
struct Data {
	int W;				// system width (grid points*)
	int H;				// system height (*)
	ptrdiff_t lW;		// width of local chunk
	ptrdiff_t lH;		// height ...
	ptrdiff_t li0;		// local x start index
	ptrdiff_t lj0;		// ... y ...
	int id;				// rank of current process
	int ID;				// number of processes
	int fold;			// order of rotational symmetry
	double dx;			// x discretization
	double dy;			// y ...
	double sk;			// kernel minor radius
	double sb;			// spread of Gaussian smoothing
	double power;		// tunable exponent for deformation field
	double *p;			// data arrays
	fftw_complex *z;
	fftw_complex *zx;
	fftw_complex *zy;
	fftw_plan z_Z;		// FFTW plans
	fftw_plan Z_z;
	fftw_plan zx_Zx;
	fftw_plan Zx_zx;
	fftw_plan Zy_zy;
	char name_n[BUFSIZ];	// density file name
	char name_phi[BUFSIZ];	// orientation file name
	char name_chi[BUFSIZ];	// deformation file name
};

// configures arrays and FFTW plans
void configure_arrays(struct Data *data) {
	int W = data->W;
	int H = data->H;
	
	ptrdiff_t A = fftw_mpi_local_size_2d_transposed(H, W, MPI_COMM_WORLD, &data->lH, &data->lj0, &data->lW, &data->li0);

	data->p = fftw_alloc_real(A);
	data->z = fftw_alloc_complex(A);
	data->zx = fftw_alloc_complex(A);
	data->zy = fftw_alloc_complex(A);
	
	fftw_complex *z = data->z;
	fftw_complex *zx = data->zx;
	fftw_complex *zy = data->zy;
	data->z_Z = fftw_mpi_plan_dft_2d(H, W, z, z, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_OUT);
	data->Z_z = fftw_mpi_plan_dft_2d(H, W, z, z, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);
	data->Zx_zx = fftw_mpi_plan_dft_2d(H, W, zx, zx, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);
	data->Zy_zy = fftw_mpi_plan_dft_2d(H, W, zy, zy, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);
}

// reads the density field
void read_state(struct Data *data) {
	char input[BUFSIZ];
	sprintf(input, "%s", data->name_n);
	FILE *file = fopen(input, "r");
	if(file == NULL) {
		printf("Error: Input file not found.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if(fscanf(file, " %d %d %lf %lf ", &data->W, &data->H, &data->dx, &data->dy) != 4) {
		printf("Error: Invalid data.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	configure_arrays(data);
	int W = data->W;
	int H = data->H;
	int lH = data->lH;
	int lj0 = data->lj0;
	double *p = data->p;
	int i, j, k;
	double dummy;
	for(j = 0; j < H; j++) {
		if(j >= lj0 && j < lj0 + lH) {
			k = W*(j - lj0);
			for(i = 0; i < W; i++) {
				if(fscanf(file, "%lf", &p[k + i]));
			}
		}
		else for(i = 0; i < W; i++) {
			if(fscanf(file, "%lf", &dummy));
		}
	}
	fclose(file);
}

// computes and writes the orientation field
void write_phi(struct Data *data) {
	int W = data->W;
	int H = data->H;
	int lW = data->lW;
	int lH = data->lH;
	int li0 = data->li0;
	int ID = data->ID;
	int fold = data->fold;
	double sk = data->sk;
	double sb = data->sb;
	double *p = data->p;
	fftw_complex *z = data->z;
	int i, gi, j, k, id;
	double dkx = 2.0*pi/(data->dx*W);
	double dky = 2.0*pi/(data->dy*H);
	double divA = 1.0/W/H;
	
	// determine global minimum
	double min = 1.0e100;
	for(j = 0; j < lH; j++) {
		k = W*j;
		for(i = 0; i < W; i++) {
			if(p[k + i] < min) min = p[k + i];
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
	// shift so that min = 0
	for(j = 0; j < lH; j++) {
		k = W*j;
		for(i = 0; i < W; i++) {
			z[k + i] = p[k + i] - min;
		}
	}
	
	double phi, kx, ky, kx2, q;
	double a = - 0.5/(sk*sk);
	
	// convolute with kernel
	fftw_execute(data->z_Z);
	for(i = 0; i < lW; i++) {
		gi = li0 + i;
		if(gi < W/2) kx = gi*dkx;
		else kx = (gi - W)*dkx;
		kx2 = kx*kx;
		k = H*i;
		for(j = 0; j < H; j++) {
			if(j < H/2) ky = j*dky;
			else ky = (j - H)*dky;
			phi = atan2(ky, kx);
			q = sqrt(kx2 + ky*ky) - 1.0;
			z[k + j] *= exp(a*q*q)*cexp(I*fold*phi)*divA;
		}
	}
	fftw_execute(data->Z_z);
	
	// multiply by shifted
	for(j = 0; j < lH; j++) {
		k = W*j;
		for(i = 0; i < W; i++) {
			z[k + i] *= p[k + i] - min;
		}
	}
	
	// smooth z
	fftw_execute(data->z_Z);
	a = - 0.5/(sb*sb);
	for(i = 0; i < lW; i++) {
		gi = li0 + i;
		if(gi < W/2) kx = gi*dkx;
		else kx = (gi - W)*dkx;
		kx2 = kx*kx;
		k = H*i;
		for(j = 0; j < H; j++) {
			if(j < H/2) ky = j*dky;
			else ky = (j - H)*dky;
			z[k + j] *= exp(a*(kx2 + ky*ky))*divA;	// Gaussian blur
			data->zx[k + j] = z[k + j];
		}
	}
	fftw_execute(data->Z_z);
	
	// normalize
	double max = - 1.0e100;
	double abs;
	for(j = 0; j < lH; j++) {
		k = W*j;
		for(i = 0; i < W; i++) {
			abs = cabs(z[k + i]);
			if(abs > max) max = abs;
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	double divmax = 1.0/max;
	for(j = 0; j < lH; j++) {
		k = W*j;
		for(i = 0; i < W; i++) {
			z[k + i] *= divmax;
		}
	}
	
	// write phi
	FILE *file;
	for(id = 0; id < ID; id++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if(id == data->id) {
			if(id == 0) {
				file = fopen(data->name_phi, "w");
				fprintf(file, "%d %d %lf %lf\n", W, H, data->dx, data->dy);
			}
			else {
				file = fopen(data->name_phi, "a");
			}
			for(j = 0; j < lH; j++) {
				k = W*j;
				for(i = 0; i < W; i++) {
					fprintf(file, "%e %e\n", creal(z[k+i]), cimag(z[k+i]));
				}
			}
			fclose(file);
		}
	}
}

// computes and writes the deformation field
void write_chi(struct Data *data) {
	int W = data->W;
	int H = data->H;
	int lW = data->lW;
	int lH = data->lH;
	int li0 = data->li0;
	int ID = data->ID;
	double dx = data->dx;
	double dy = data->dy;
	double *p = data->p;
	fftw_complex *z = data->z;
	fftw_complex *zx = data->zx;
	fftw_complex *zy = data->zy;
	int i, gi, j, k, id;
	double kx, kx2, ky;
	double dkx = 2.0*pi/(dx*W);
	double dky = 2.0*pi/(dy*H);
	double divA = 1.0/W/H;
	
	// clear p
	for(j = 0; j < lH; j++) {
		k = W*j;
		for(i = 0; i < W; i++) {
			p[k + i] = 0.0;
		}
	}
	
	// restore Z
	for(i = 0; i < lW; i++) {
		k = H*i;
		for(j = 0; j < H; j++) {
			z[k + j] = zx[k + j];
		}
	}
	
	// compute dphi/dx and dphi/dy
	for(i = 0; i < lW; i++) {
		gi = li0 + i;
		if(gi < W/2) kx = gi*dkx;
		else kx = (gi - W)*dkx;
		k = H*i;
		for(j = 0; j < H; j++) {
			if(j < H/2) ky = j*dky;
			else ky = (j - H)*dky;
			zx[k + j] = I*kx*z[k + j];
			zy[k + j] = I*ky*z[k + j];
		}
	}
	fftw_execute(data->Zx_zx);
	fftw_execute(data->Zy_zy);
	
	// compute |grad(phi)|^p
	double u, v;
	for(j = 0; j < lH; j++) {
		k = W*j;
		for(i = 0; i < W; i++) {
			u = cabs(zx[k + i]);
			v = cabs(zy[k + i]);
			z[k + i] = pow(u*u + v*v, 0.5*data->power);
		}
	}
	
	// compute the sum of normalized convolutions
	fftw_execute(data->z_Z);
	
	double l0 = 4.0*pi/sqrt(3.0);
	double L = W*dx;
	if(H*dy < L) L = H*dy;
	int n;
	for(n = 2; l0*n < 0.5*L; n *= 2) {
		
		// convolve |grad(phi)|^p with Gaussian kernel
		double s = 1.0/n;
		double a = - 0.5/(s*s);
		for(i = 0; i < lW; i++) {
			gi = li0 + i;
			if(gi < W/2) kx = gi*dkx;
			else kx = (gi - W)*dkx;
			kx2 = kx*kx;
			k = H*i;
			for(j = 0; j < H; j++) {
				if(j < H/2) ky = j*dky;
				else ky = (j - H)*dky;
				zx[k + j] = exp(a*(kx2 + ky*ky))*divA*z[k + j];
			}
		}
		fftw_execute(data->Zx_zx);
		
		// determine global maximum
		double max = - 1.0e100;
		for(j = 0; j < lH; j++) {
			k = W*j;
			for(i = 0; i < W; i++) {
				if(creal(zx[k + i]) > max) max = creal(zx[k + i]);
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		
		// normalize
		double divmax = 1.0/max;
		for(j = 0; j < lH; j++) {
			k = W*j;
			for(i = 0; i < W; i++) {
				p[k + i] += creal(zx[k + i])*divmax;
			}
		}
		
	}
	
	// write chi
	FILE *file;
	for(id = 0; id < ID; id++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if(id == data->id) {
			if(id == 0) {
				file = fopen(data->name_chi, "w");
				fprintf(file, "%d %d %lf %lf\n", W, H, data->dx, data->dy);
			}
			else {
				file = fopen(data->name_chi, "a");
			}
			for(j = 0; j < lH; j++) {
				k = W*j;
				for(i = 0; i < W; i++) {
					fprintf(file, "%e\n", p[k + i]);
				}
			}
			fclose(file);
		}
	}
}

// frees allocated memory
void clear_arrays(struct Data *data) {
	fftw_free(data->p);
	fftw_free(data->zx);
	fftw_free(data->zy);
}

int main(int argc, char **argv) {
	// init MPI
	MPI_Init(&argc, &argv);
	fftw_mpi_init();
	
	// create structs
	struct Data data;
	MPI_Comm_rank(MPI_COMM_WORLD, &data.id);
	if(argc == 1) {			// print help
		if(data.id == 0) {
			printf("This is a tool for computing the orientation and deformation fields corresponding to a density field describing a periodic system such as a polycrystalline microstructure. For more information, see Hirvonen et al., Physical Review Materials 2, 103603 (2018), DOI: 10.1103/PhysRevMaterials.2.103603 or arXiv:1806.00700.\n\n");
			
			printf("Examples of valid syntax:\n\n");
			
			printf("mpirun -np N fields density.n orientation.phi [order of rotational symmetry] [spread of complex-valued kernel] [spread of Gaussian smoothing kernel]\n\n");
			
			printf("mpirun -np N fields density.n orientation.phi deformation.chi [order of rotational symmetry] [spread of complex-valued kernel] [spread of Gaussian smoothing kernel] [tunable exponent for deformation field]\n\n");
			
			printf("The former example computes only the orientation field, whereas the latter computes also the deformation field. The program expects the names of the input and output files, the order of rotational symmetry (2 for stripe, 4 for square, 6 for hexagonal or honeycomb, and n for n-fold quasilattices; note that the method works only for even-fold symmetries), spread of the complex-valued kernel (in reciprocal units*), spread of the Gaussian smoothing kernel (*) and the tunable exponent for the deformation field. Furthermore, the input field is expected to have an approximate characteristic length scale k = 1 in reciprocal units.\n\n");
		}
		MPI_Finalize();
		return 0;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &data.ID);
	
	if(argc != 6 && argc != 8) {
		if(data.id == 0) {
			printf("Error: Invalid syntax.\n");
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	if(argc == 6) {
		strcpy(data.name_n, argv[1]);	// density file name
		strcpy(data.name_phi, argv[2]);	// orientation file name
		data.fold = atoi(argv[3]);		// order of rotational symmetry
		data.sk = atof(argv[4]);		// kernel minor radius
		data.sb = atof(argv[5]);		// spread of Gaussian smoothing
	}
	else if(argc == 8) {
		strcpy(data.name_n, argv[1]);	// density file name
		strcpy(data.name_phi, argv[2]);	// orientation file name
		strcpy(data.name_chi, argv[3]);	// deformation file name
		data.fold = atoi(argv[4]);		// order of rotational symmetry
		data.sk = atof(argv[5]);		// kernel minor radius
		data.sb = atof(argv[6]);		// spread of Gaussian smoothing
		data.power = atof(argv[7]);		// tunable exponent for deformation field
	}
	
	read_state(&data);
	write_phi(&data);
	if(argc == 8) write_chi(&data);
	
	clear_arrays(&data);
	MPI_Finalize();

	return 0;

}
