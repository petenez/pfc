// Author: Petri Hirvonen, petenez@gmail.com, 5 September 2019

// compile: mpicc pfc-relax.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o pfc-relax

// print instructions: mpirun -np 1 pfc-relax

// This program implements a multimode and multicomponent PFC model. It is given by a free energy

// F = \int d \boldsymbol{r} \left( \sum_{u = 0}^{U - 1} \left( \lambda_u n_u + \frac{\alpha_u}{2} n_u^2 + \frac{\beta_u}{2} n_u \left( \prod_{w = 0}^{W - 1} \left( \omega_{uw} \left( \nu_{uw}^2 + \nabla^2 \right)^2 + \varphi_{uw} \right) \right) n_u + \frac{\gamma_u}{3} n_u^3 + \frac{\delta_u}{4} n_u^4 \right) + \sum_{u = 0}^{U - 2} \sum_{v = u + 1}^{U - 1} \left( \alpha_{uv} n_u n_v + \beta_{uv} n_u \left( \mu_{uv}^2 + \nabla^2 \right)^2 n_v + \frac{\gamma_{uv}}{2} \left( n_u^2 n_v + n_u n_v^2 \right) + \varepsilon_{uv} \eta_u \eta_v \right) \right),

// where n_u is the uth density field, \eta_u = G \ast n_u where G is a Gaussian smoothing kernel with the Fourier transform \hat{G} \left( \boldsymbol{k} \right) = e^{-\left| \boldsymbol{k} \right|^2 / \left( 2 \sigma^2 \right)} and the asterisk denotes a convolution. The dynamics of the density fields are given by

// \frac{\partial n_u}{\partial t} = \square \frac{\delta F}{\delta n_u} = \square \left( \lambda_u + \alpha_u n_u + \beta_u \left( \prod_{w = 0}^{W - 1} \left( \omega_{uw} \left( \nu_{uw}^2 + \nabla^2 \right)^2 + \varphi_{uw} \right) \right) n_u + \gamma_u n_u^2 + \delta_u n_u^3 + \sum_{v = 0, v \neq u}^{U - 1} \left( \alpha_{uv} n_v + \beta_{uv} \left( \mu_{uv}^2 + \nabla^2 \right)^2 n_v + \frac{\gamma_{uv}}{2} \left( 2 n_u n_v + n_v^2 \right) + \varepsilon_{uv} G \ast \eta_v \right) \right),

// where the \square = \nabla^2 for conserved dynamics and \square = -1 for nonconserved dynamics. The functional derivative of F with respect to n_u is given by \frac{\delta F}{\delta n_u}. The dynamics can be written as

// \frac{\partial n_u}{\partial t} = A n_u + B m_u,

// where A is a "linear operator" acting on n_u and B is a "nonlinear operator" acting on the nonlinear parts. The dynamics can be solved numerically using the semi-implicit spectral method given in the book Phase-field methods in material science and engineering by Provatas and Elder (http://www.physics.mcgill.ca/~provatas/papers/Phase_Field_Methods_text.pdf).

// This code follows the above notation. For more information about PFC, see my thesis (available: http://urn.fi/URN:ISBN:978-952-60-8608-8) and our recent paper on PFC heterostructures (arXiv:1908.05564), and the references therein.

#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define pi 3.14159265358979323846264338327


struct S {

	// ranks
	int r;					// a process's rank
	int R;					// # of ranks
	
	int t;					// relaxation step
	int t0;					// relaxation start wall-clock time
	int tend;				// relaxation last step
	
	// calculation grid dimensions etc
	int Nx;					// grid points in x direction
	int Ny;					// grid points in y direction
	
	// local chunk dimensions and start indices
	ptrdiff_t lNx;			// grid points in x direction in local chunk
	ptrdiff_t lNy;			// grid points in y direction in local chunk
	ptrdiff_t li0;			// local x start index
	ptrdiff_t lj0;			// local y start index
	
	// input and output
	int Tprint;				// interval for printing output
	int Twrite;				// interval for writing out data
	int pprint;				// precision for printing output
	int pwrite;				// precision for writing out data
	int* twrite;			// relaxation steps for writing out data
	
	char input[BUFSIZ];		// input filename
	char* fname;			// free energy density output filename
	
	char** innames;			// density input filenames
	char** outnames;		// density output filenames
	
	// discretization
	double dx;				// x discretization
	double dy;				// y discretization
	double dt;				// time step size
	
	// system size optimization
	int Topt;				// interval for system size optimization
	int xopt;				// system size optimized in x direction?
	int yopt;				// system size optimized in y direction?
	double d;				// system size optimization step size
	
	// parameter arrays
	int U;					// # of density fields
	int* W;					// # of modes for each density field
	double c;				// dynamics (0: nonconserved, !0: conserved)
	double sigma;			// see the model definition above
	
	double* lambda;			// see ...
	double* delta;			// ...
	
	double** alpha;
	double** beta;
	double** gamma;
	double** epsilon;
	double** nu2;
	double** omega;
	double** phi;
	double** mu2;
	
	// data arrays
	double* mc;				// array for nonlinear couplings
	
	double fave;			// average free energy density
	
	double* nave;			// average densities
	
	double* f;				// free energy density
	
	double** A;				// arrays for operators A
	double** B;				// arrays for operators B
	double** n;				// arrays for densities
	double** m;				// arrays for nonlinear parts of densities
	double** e;				// arrays for smoothed densities
	double*** q;			// arrays for gradient term couplings
	
	// FFTW plans
	// lowercase: direct space; uppercase: Fourier space
	fftw_plan* n_N;
	fftw_plan* m_M;
	fftw_plan* M_m;
	fftw_plan* E_e;
	fftw_plan** Q_q;

};

// allocates parameter arrays
void allocate_parameters(struct S* s) {

	int u;
	int U = s->U;

	// W, lambda, delta
	s->W = fftw_malloc(U*sizeof(int));
	s->lambda = fftw_malloc(U*sizeof(double));
	s->delta = fftw_malloc(U*sizeof(double));
	
	// nu2, omega, phi
	s->nu2 = fftw_malloc(U*sizeof(double*));
	s->omega = fftw_malloc(U*sizeof(double*));
	s->phi = fftw_malloc(U*sizeof(double*));
	// second dimensions unknown --> left uninitialized for now
	
	// alpha, beta, gamma, epsilon, mu2
	s->alpha = fftw_malloc(U*sizeof(double*));
	s->beta = fftw_malloc(U*sizeof(double*));
	s->gamma = fftw_malloc(U*sizeof(double*));
	s->epsilon = fftw_malloc(U*sizeof(double*));
	s->mu2 = fftw_malloc(U*sizeof(double*));
	int U2 = U*U;
	s->alpha[0] = fftw_malloc(U2*sizeof(double));
	s->beta[0] = fftw_malloc(U2*sizeof(double));
	s->gamma[0] = fftw_malloc(U2*sizeof(double));
	s->epsilon[0] = fftw_malloc(U2*sizeof(double));
	s->mu2[0] = fftw_malloc(U2*sizeof(double));
	for(u = 1; u < U; u++) {
		s->alpha[u] = s->alpha[u - 1] + U;
		s->beta[u] = s->beta[u - 1] + U;
		s->gamma[u] = s->gamma[u - 1] + U;
		s->epsilon[u] = s->epsilon[u - 1] + U;
		s->mu2[u] = s->mu2[u - 1] + U;
	}
	
}

// allocates data arrays
void allocate_data(struct S* s) {

	int u, v;
	int U = s->U;
	int lNxNyh = fftw_mpi_local_size_2d_transposed(s->Ny, s->Nx/2+1, MPI_COMM_WORLD, &s->lNy, &s->lj0, &s->lNx, &s->li0);
	int lNxNyp = 2*lNxNyh;

	// mc, nave
	s->mc = fftw_malloc(U*sizeof(double));
	s->nave = fftw_malloc(U*sizeof(double));
	
	// f, A, B, n, m, e, q
	s->f = fftw_malloc(lNxNyp*sizeof(double));
	s->A = fftw_malloc(U*sizeof(double*));
	s->B = fftw_malloc(U*sizeof(double*));
	s->n = fftw_malloc(U*sizeof(double*));
	s->m = fftw_malloc(U*sizeof(double*));
	s->e = fftw_malloc(U*sizeof(double*));
	s->q = fftw_malloc(U*sizeof(double**));
	for(u = 0; u < U; u++) {
		s->A[u] = fftw_malloc(lNxNyh*sizeof(double));
		s->B[u] = fftw_malloc(lNxNyh*sizeof(double));
		s->n[u] = fftw_malloc(lNxNyp*sizeof(double));
		s->m[u] = fftw_malloc(lNxNyp*sizeof(double));
		s->e[u] = fftw_malloc(lNxNyp*sizeof(double));
		s->q[u] = fftw_malloc(U*sizeof(double*));
		for(v = u + 1; v < U; v++) {
			s->q[u][v] = fftw_malloc(lNxNyp*sizeof(double));
		}
	}

}

// sets up FFTW plans (real-data transforms with transposed order in Fourier space)
void plan_FFTs(struct S* s) {

	int u, v;
	int U = s->U;
	int Nx = s->Nx;
	int Ny = s->Ny;
	
	double** n = s->n;
	double** m = s->m;
	double** e = s->e;
	double*** q = s->q;
	
	fftw_complex** N = (fftw_complex**)n;
	fftw_complex** M = (fftw_complex**)m;
	fftw_complex** E = (fftw_complex**)e;
	fftw_complex*** Q = (fftw_complex***)q;

	s->n_N = fftw_malloc(U*sizeof(fftw_plan));
	s->m_M = fftw_malloc(U*sizeof(fftw_plan));
	s->M_m = fftw_malloc(U*sizeof(fftw_plan));
	s->E_e = fftw_malloc(U*sizeof(fftw_plan));
	s->Q_q = fftw_malloc(U*sizeof(fftw_plan*));
	
	fftw_plan* n_N = s->n_N;
	fftw_plan* m_M = s->m_M;
	fftw_plan* M_m = s->M_m;
	fftw_plan* E_e = s->E_e;
	fftw_plan** Q_q = s->Q_q;
	
	for(u = 0; u < U; u++) {
		// FFT
		n_N[u] = fftw_mpi_plan_dft_r2c_2d(Ny, Nx, n[u], N[u], MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);
		m_M[u] = fftw_mpi_plan_dft_r2c_2d(Ny, Nx, m[u], M[u], MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);
		
		// FFT^(-1)
		M_m[u] = fftw_mpi_plan_dft_c2r_2d(Ny, Nx, M[u], m[u], MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_IN);
		E_e[u] = fftw_mpi_plan_dft_c2r_2d(Ny, Nx, E[u], e[u], MPI_COMM_WORLD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);
		Q_q[u] = fftw_malloc(U*sizeof(fftw_plan));
		for(v = u + 1; v < U; v++) {
			Q_q[u][v] = fftw_mpi_plan_dft_c2r_2d(Ny, Nx, Q[u][v], q[u][v], MPI_COMM_WORLD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);
		}
	}

}

// prints an error message and terminates the program
void fail(struct S* s, int i) {

	switch(i) {
	
		case 0:
			if(s->r == 0) printf("Error: Input file not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 1:
			if(s->r == 0) printf("Error: Invalid multicharacter label found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
		
		case 2:
			if(s->r == 0) printf("Error: Failed to read an input data file name.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 3:
			if(s->r == 0) printf("Error: Input data file not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 4:
			if(s->r == 0) printf("Error: Invalid input data file header.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 5:
			if(s->r == 0) printf("Error: Invalid mismatching dimensions of input data files.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 6:
			if(s->r == 0) printf("Error: Invalid number of output data files specified.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 7:
			if(s->r == 0) printf("Error: Failed to read an output data file name.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 8:
			if(s->r == 0) printf("Error: Failed to read the free energy density output data file name.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 9:
			if(s->r == 0) printf("Error: Invalid parameters for writing output.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 10:
			if(s->r == 0) printf("Error: Invalid write steps specified.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 11:
			if(s->r == 0) printf("Error: Failed to read the number of density fields.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 12:
			if(s->r == 0) printf("Error: Mismatching numbers of input data files and density fields specified.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 13:
			if(s->r == 0) printf("Error: Label for parameters lambda not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 14:
			if(s->r == 0) printf("Error: Invalid syntax for parameters lambda.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 15:
			if(s->r == 0) printf("Error:  Label for parameters alpha not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 16:
			if(s->r == 0) printf("Error: Invalid syntax for parameters alpha.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 17:
			if(s->r == 0) printf("Error:  Label for parameters beta not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 18:
			if(s->r == 0) printf("Error: Invalid syntax for parameters beta.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 19:
			if(s->r == 0) printf("Error:  Label for parameters gamma not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 20:
			if(s->r == 0) printf("Error: Invalid syntax for parameters gamma.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 21:
			if(s->r == 0) printf("Error: Label for parameters delta not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 22:
			if(s->r == 0) printf("Error: Invalid syntax for parameters delta.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 23:
			if(s->r == 0) printf("Error: Label for parameters epsilon not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 24:
			if(s->r == 0) printf("Error: Invalid syntax for parameters epsilon.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 25:
			if(s->r == 0) printf("Error: Label for parameters sigma not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 26:
			if(s->r == 0) printf("Error: Invalid syntax for parameters sigma.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 27:
			if(s->r == 0) printf("Error: Label for parameters 1/nu not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 28:
			if(s->r == 0) printf("Error: Invalid syntax for parameters 1/nu.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 29:
			if(s->r == 0) printf("Error: Label for parameters omega not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 30:
			if(s->r == 0) printf("Error: Invalid syntax for parameters omega.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 31:
			if(s->r == 0) printf("Error: Label for parameters phi not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 32:
			if(s->r == 0) printf("Error: Invalid syntax for parameters phi.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 33:
			if(s->r == 0) printf("Error: Label for parameters 1/mu not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 34:
			if(s->r == 0) printf("Error: Invalid syntax for parameters 1/mu.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 35:
			if(s->r == 0) printf("Error: Invalid parameters for relaxation.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 36:
			if(s->r == 0) printf("Error: Invalid label found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 37:
			if(s->r == 0) printf("Error: Input data file not found.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 38:
			if(s->r == 0) printf("Error: Invalid data in data file.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 39:
			if(s->r == 0) printf("Error: Not a number value detected.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
	
		case 40:
			if(s->r == 0) printf("Error: An error occurred while reading the input file.\n");
			MPI_Abort(MPI_COMM_WORLD, i);
		
	}
}

// reads the input file
void read_input(struct S* s) {

	FILE* input = fopen(s->input, "r");
	if(input == NULL) fail(s, 0);

	int u, v, w, U, V, Nx, Ny, Nx_, Ny_, N, M;
	U = 0;
	int* W;
	int offset;
	double tmp;
	double** nu2;
	double** mu2;
	char str[BUFSIZ];
	char dummy[BUFSIZ];
	char* fname = NULL;
	char* pstr;
	
	while(fscanf(input, " %s", str) > 0) {

		// comment
		if(str[0] == '#') {
			if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
		}

		// garbage
		else if(strlen(str) > 1) fail(s, 1);

		// input and output data files and write parameters
		else if(str[0] == 'I') {

			// input data files
			if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
			pstr = str;
			while(sscanf(pstr, "%s%n", dummy, &offset) == 1) {
				U++;
				pstr += offset;
			}
			s->U = U;	// # of input files
			char** innames = fftw_malloc(U*sizeof(char*));
			for(u = 0; u < U; u++) innames[u] = fftw_malloc(BUFSIZ*sizeof(char));
			s->innames = innames;
			pstr = str;
			Nx = -1;
			Ny = -1;
			// make sure dimensions match
			for(u = 0; u < U; u++) {
				if(sscanf(pstr, "%s%n", innames[u], &offset) != 1) fail(s, 2);
				pstr += offset;
				
				FILE* file = fopen(innames[u], "r");
				if(file == NULL) fail(s, 3);
				if(fscanf(file, "%d %d %lf %lf", &Nx_, &Ny_, &tmp, &tmp) != 4) fail(s, 4);
				if((Nx > 0 && Nx != Nx_) || (Ny > 0 && Ny != Ny_)) fail(s, 5);
				Nx = Nx_;
				Ny = Ny_;
				fclose(file);
			}
			s->Nx = Nx;
			s->Ny = Ny;

			// output data files
			if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
			pstr = str;
			V = 0;
			// make sure # of output files matches 
			while(sscanf(pstr, "%s%n", dummy, &offset) == 1) {
				V++;
				pstr += offset;
			}
			if(V < U || V > U + 1) fail(s, 6);
			char** outnames = fftw_malloc(U*sizeof(char*));
			for(u = 0; u < U; u++) outnames[u] = fftw_malloc(BUFSIZ*sizeof(char));
			s->outnames = outnames;
			pstr = str;
			// density
			for(u = 0; u < U; u++) {
				if(sscanf(pstr, "%s%n", outnames[u], &offset) != 1) fail(s, 7);
				pstr += offset;
			}
			// free energy density
			if(V > U) {
				fname = fftw_malloc(BUFSIZ*sizeof(char));
				s->fname = fname;
				if(sscanf(pstr, "%s", fname) != 1) fail(s, 8);
			}
			s->fname = fname;

			// write parameters
			s->Twrite = -1;
			if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
			pstr = str;
			N = 0;
			while(sscanf(pstr, "%s%n", dummy, &offset) == 1) {
				N++;
				pstr += offset;
			}
			// write interval specified
			if(N == 4) {
				if(sscanf(str, " %d %d %d %d", &s->Tprint, &s->Twrite, &s->pprint, &s->pwrite) != 4) fail(s, 9);
			}
			// write steps specified
			else if(N >= 6) {
				pstr = str;
				M = N - 5;
				int* twrite = fftw_malloc(M*sizeof(int));
				s->twrite = twrite;
				if(sscanf(pstr, "%d%n", &s->Tprint, &offset) != 1) fail(s, 9);
				pstr += offset;
				if(sscanf(pstr, "%s%n", dummy, &offset) != 1) fail(s, 9);
				pstr += offset;
				if(strcmp(dummy, "[") != 0) fail(s, 9);
				for(u = 0; u < M; u++) {
					if(sscanf(pstr, "%d%n", &twrite[u], &offset) != 1) fail(s, 9);
					pstr += offset;
				}
				for(u = 1; u < M; u++) {
					if(twrite[u] <= twrite[u - 1]) fail(s, 10);
				}
				if(sscanf(pstr, "%s%n", dummy, &offset) != 1) fail(s, 9);
				pstr += offset;
				if(strcmp(dummy, "]") != 0) fail(s, 9);
				if(sscanf(pstr, "%d %d", &s->pprint, &s->pwrite) != 2) fail(s, 9);
			}
			else fail(s, 9);

		}

		// model parameters
		else if(str[0] == 'M') {

			// # of density fields
			if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
			if(sscanf(str, "%d", &V) != 1) fail(s, 11);
			if(U != V) fail(s, 12);
			allocate_parameters(s);

			// lambda
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 13);
				if(strcmp(str, "l") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				if(sscanf(str, " %lf", &s->lambda[u]) != 1) fail(s, 14);
			}

			// alpha
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 15);
				if(strcmp(str, "a") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				for(v = u; v < U; v++) {
					if(sscanf(pstr, " %lf%n", &s->alpha[u][v], &offset) != 1) fail(s, 16);
					pstr += offset;
				}
			}

			// beta
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 17);
				if(strcmp(str, "b") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				for(v = u; v < U; v++) {
					if(sscanf(pstr, " %lf%n", &s->beta[u][v], &offset) != 1) fail(s, 18);
					pstr += offset;
				}
			}

			// gamma
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 19);
				if(strcmp(str, "g") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				for(v = u; v < U; v++) {
					if(sscanf(pstr, " %lf%n", &s->gamma[u][v], &offset) != 1) fail(s, 20);
					pstr += offset;
				}
			}

			// delta
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 21);
				if(strcmp(str, "d") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				if(sscanf(str, " %lf", &s->delta[u]) != 1) fail(s, 22);
			}

			// epsilon
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 23);
				if(strcmp(str, "e") == 0) break;
			}
			for(u = 0; u < U - 1; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				for(v = u + 1; v < U; v++) {
					if(sscanf(pstr, " %lf%n", &s->epsilon[u][v], &offset) != 1) fail(s, 24);
					pstr += offset;
				}
			}

			// sigma
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 25);
				if(strcmp(str, "s") == 0) break;
			}
			if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
			if(sscanf(str, " %lf", &s->sigma) != 1) fail(s, 26);

			// 1/nu
			nu2 = s->nu2;
			W = s->W;
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 27);
				if(strcmp(str, "n") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				W[u] = 0;
				while(sscanf(pstr, "%s%n", dummy, &offset) == 1) {
					W[u]++;
					pstr += offset;
				}
				nu2[u] = fftw_malloc(W[u]*sizeof(double));
				s->omega[u] = fftw_malloc(W[u]*sizeof(double));
				s->phi[u] = fftw_malloc(W[u]*sizeof(double));
				pstr = str;
				for(w = 0; w < W[u]; w++) {
					if(sscanf(pstr, "%lf%n", &nu2[u][w], &offset) != 1) fail(s, 28);
					nu2[u][w] = 1.0/nu2[u][w]/nu2[u][w];
					pstr += offset;
				}
			}

			// omega
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 29);
				if(strcmp(str, "o") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				for(w = 0; w < W[u]; w++) {
					if(sscanf(pstr, "%lf%n", &s->omega[u][w], &offset) != 1) fail(s, 30);
					pstr += offset;
				}
			}

			// phi
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 31);
				if(strcmp(str, "p") == 0) break;
			}
			for(u = 0; u < U; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				for(w = 0; w < W[u]; w++) {
					if(sscanf(pstr, "%lf%n", &s->phi[u][w], &offset) != 1) fail(s, 32);
					pstr += offset;
				}
			}

			// 1/mu
			mu2 = s->mu2;
			while(1) {
				if(fscanf(input, " %s", str) != 1) fail(s, 33);
				if(strcmp(str, "m") == 0) break;
			}
			for(u = 0; u < U - 1; u++) {
				if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
				pstr = str;
				for(v = u + 1; v < U; v++) {
					if(sscanf(pstr, " %lf%n", &mu2[u][v], &offset) != 1) fail(s, 34);
					if(s->beta[u][v] == 0.0) mu2[u][v] = 0.0;
					else mu2[u][v] = 1.0/mu2[u][v]/mu2[u][v];
					pstr += offset;
				}
			}
			
		}

		// relaxation
		else if(str[0] == 'R') {
			if(fgets(str, BUFSIZ, input) == NULL) fail(s, 40);
			if(sscanf(str, "%lf %d %lf %lf %lf %d %d %d", &s->c, &s->tend, &s->dx, &s->dy, &s->dt, &s->Topt, &s->xopt, &s->yopt) != 8) fail(s, 35);
			
			if(s->Twrite > 0) {
				N = s->tend/s->Twrite + 1;
				s->twrite = fftw_malloc(N*sizeof(int));
				for(u = 0; u < N; u++) {
					s->twrite[u] = s->Twrite*u;
				}
			}
		}

		// garbage
		else fail(s, 36);
	}

	fclose(input);
}

// read input data files to initialize density fields
void read_data(struct S* s) {
	int U = s->U;
	int Nx = s->Nx;
	int Ny = s->Ny;
	int Nxp = 2*(Nx/2 + 1);
	int lNy = s->lNy;
	int lj0 = s->lj0;
	char** innames = s->innames;
	double** m = s->m;
	int i, j, ij, u;
	int iummy;					// dummy variables
	double dummy;
	FILE* file;
	// different density field
	for(u = 0; u < U; u++) {
		file = fopen(innames[u], "r");
		if(file == NULL) fail(s, 37);	// file not found
		// read dimensions and discretizations
		if(fscanf(file, " %d %d %lf %lf ", &iummy, &iummy, &dummy, &dummy) != 4) fail(s, 4);
		for(j = 0; j < Ny; j++) {
			if(j >= lj0 && j < lj0 + lNy) {	// withing local chunk?
				ij = Nxp*(j - lj0);			// start index of current row
				for(i = 0; i < Nx; i++) {
					if(fscanf(file, " %lf ", &m[u][ij]) != 1) fail(s, 38);
					ij++;
				}
			}
			else {	// else flush from input
				for(i = 0; i < Nx; i++) {
					if(fscanf(file, " %lf ", &dummy) != 1) fail(s, 38);
				}
			}
		}
		fclose(file);
	}
}

// writes out density fields
void write_data(struct S* s) {
	int U = s->U;
	int r = s->r;
	int R = s->R;
	int Nx = s->Nx;
	int Ny = s->Ny;
	int Nxp = 2*(Nx/2 + 1);
	int lNy = s->lNy;
	int pwrite = s->pwrite;
	double dx = s->dx;
	double dy = s->dy;
	int t = s->t;
	char** outnames = s->outnames;
	double** m = s->m;
	int i, j, ij, u, p;
	char* x;
	char str[BUFSIZ];
	FILE* file;
	// different density fields
	for(u = 0; u < U; u++) {
		memset(str, '\0', BUFSIZ);
		x = strchr(outnames[u], '#');
		if(x == NULL) strcpy(str, outnames[u]);
		else {	// replace "#" with current relaxation step
			strncpy(str, outnames[u], x - outnames[u]);
			sprintf(str, "%s%d%s", str, t, x + 1);
		}
		// different processes
		for(p = 0; p < R; p++) {
			MPI_Barrier(MPI_COMM_WORLD);
			if(r == p) {
				if(r == 0) {
					file = fopen(str, "w");		// first process overwrites
					fprintf(file, "%d %d %.*lf %.*lf\n", Nx, Ny, pwrite, dx, pwrite, dy);
				}
				else file = fopen(str, "a");	// others append
				for(j = 0; j < lNy; j++) {
					ij = Nxp*j;
					for(i = 0; i < Nx; i++) {
						fprintf(file, "%.*e\n", pwrite, m[u][ij]);
						ij++;
					}
				}
				fclose(file);
			}
		}
	}
}

// writes out free energy density (see previous for more details)
void write_energy(struct S* s) {
	int r = s->r;
	int R = s->R;
	int Nx = s->Nx;
	int Ny = s->Ny;
	int Nxp = 2*(Nx/2 + 1);
	int lNy = s->lNy;
	int pwrite = s->pwrite;
	double dx = s->dx;
	double dy = s->dy;
	int t = s->t;
	char* fname = s->fname;
	double* f = s->f;
	int i, j, ij, p;
	char* x;
	char str[BUFSIZ];
	FILE* file;
	memset(str, '\0', BUFSIZ);
	x = strchr(fname, '#');
	if(x == NULL) strcpy(str, fname);
	else {
		strncpy(str, fname, x - fname);
		sprintf(str, "%s%d%s", str, t, x + 1);
	}
	for(p = 0; p < R; p++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if(r == p) {
			if(r == 0) {
				file = fopen(str, "w");
				fprintf(file, "%d %d %.*lf %.*lf\n", Nx, Ny, pwrite, dx, pwrite, dy);
			}
			else file = fopen(str, "a");
			for(j = 0; j < lNy; j++) {
				ij = Nxp*j;
				for(i = 0; i < Nx; i++) {
					fprintf(file, "%.*e\n", pwrite, f[ij]);
					ij++;
				}
			}
			fclose(file);
		}
	}
}

// prints output
void print_output(struct S* s) {
	if(s->r != 0) return;	// only 0th process
	
	int u;
	int t0 = s->t0;
	int t = s->t;
	int pprint = s->pprint;
	int U = s->U;
	
	double dx = s->dx;
	double dy = s->dy;
	double fave = s->fave;
	double* nave = s->nave;
	
	char line[BUFSIZ];
	char* input = s->input;
	char str[BUFSIZ];
	FILE* output;
	
	sprintf(line, "%d %d %.*lf %.*lf %.*lf ", t, (int)(time(NULL) - t0), pprint, dx, pprint, dy, pprint, fave);
	for(u = 0; u < U; u++) {
		sprintf(line, "%s%.*lf ", line, pprint, nave[u]);
	}
	printf("%s\n", line);
	
	memset(str, '\0', BUFSIZ);
	// input file ends with ".in"? --> replace with ".out"
	if(strcmp(input + strlen(input) - 3, ".in") == 0) {
		strncpy(str, input, strlen(input) - 3);
	}
	else strcpy(str, input);	// append ".out"
	strcat(str, ".out");
	if(s->t == 0) output = fopen(str, "w");
	else output = fopen(str, "a");
	fprintf(output, "%s\n", line);
	fclose(output);
}

// updates operators A and B
void update_AB(struct S* s) {
	int i, j, ij, gi, u, w;
	
	int U = s->U;
	int* W = s->W;
	int Nx = s->Nx;
	int Ny = s->Ny;
	int lNx = s->lNx;
	int li0 = s->li0;

	double kx2, ky, k2, P, p, a, o, ex;
	double _NxNy = 1.0/Nx/Ny;
	double dx = s->dx;
	double dy = s->dy;
	double dt = s->dt;
	double dkx = 2.0*pi/(dx*Nx);
	double dky = 2.0*pi/(dy*Ny);
	
	double c = s->c;
	
	double** nu2 = s->nu2;
	double** omega = s->omega;
	double** phi = s->phi;
	double** alpha = s->alpha;
	double** beta = s->beta;
	
	double** A = s->A;
	double** B = s->B;
	
	if(c != 0.0) c = 1.0;
	for(i = 0; i < lNx; i++) {
		gi = li0 + i;
		kx2 = gi*dkx;
		kx2 *= kx2;
		ij = Ny*i;
		for(j = 0; j < Ny; j++) {
			if(j < Ny/2) ky = j*dky;
			else ky = (j - Ny)*dky;
			k2 = kx2 + ky*ky;
			for(u = 0; u < U; u++) {
				// multiple modes
				P = 1.0;
				for(w = 0; w < W[u]; w++) {
					p = nu2[u][w] - k2;
					P *= omega[u][w]*p*p + phi[u][w];
				}
				a = alpha[u][u] + beta[u][u]*P;
				o = c - 1.0 - c*k2;
				ex = exp(o*a*dt);
				A[u][ij] = ex;
				if(a == 0.0) B[u][ij] = o*dt;	// avoid divide by zero
				else B[u][ij] = (ex - 1.0)/a;
				B[u][ij] *= _NxNy;
			}
			ij++;
		}
	}
}

void initialize(struct S* s) {
	allocate_data(s);
	plan_FFTs(s);
	read_data(s);
	update_AB(s);
}

// descales FFT'd fields
void descale(struct S* s, double** a) {
	int U = s->U;
	int Nx = s->Nx;
	int Ny = s->Ny;
	int lNx = s->lNx;
	double _NxNy = 1.0/Nx/Ny;
	int i, j, ij, u;
	fftw_complex** A = (fftw_complex**)a;
	for(u = 0; u < U; u++) {
		for(i = 0; i < lNx; i++) {
			ij = Ny*i;
			for(j = 0; j < Ny; j++) {
				A[u][ij] *= _NxNy;
				ij++;
			}
		}
	}
}

// computes the free energy of the system and the average densities of its density fields
void evaluate(struct S* s) {
	
	int i, j, ij, gi, u, v, w;
	
	int U = s->U;
	int Nx = s->Nx;
	int Ny = s->Ny;
	int Nxp = 2*(Nx/2 + 1);
	int lNx = s->lNx;
	int lNy = s->lNy;
	int li0 = s->li0;
	
	int* W = s->W;
	
	double a, kx2, ky, k2, G, d, P, p, nunu, nvnv;
	double _NxNy = 1.0/Nx/Ny;
	double dx = s->dx;
	double dy = s->dy;
	double dkx = 2.0*pi/(dx*Nx);
	double dky = 2.0*pi/(dy*Ny);
	double _3 = 1.0/3.0;
	
	double sigma = s->sigma;
	double* lambda = s->lambda;
	double* delta = s->delta;
	double** alpha = s->alpha;
	double** beta = s->beta;
	double** gamma = s->gamma;
	double** epsilon = s->epsilon;
	double** nu2 = s->nu2;
	double** omega = s->omega;
	double** phi = s->phi;
	double** mu2 = s->nu2;
	
	double* fave = &s->fave;
	double* nave = s->nave;
	
	double* f = s->f;
	double** n = s->n;
	double** m = s->m;
	double** e = s->e;
	double*** q = s->q;
	
	fftw_complex** M = (fftw_complex**)m;
	fftw_complex** E = (fftw_complex**)e;
	fftw_complex*** Q = (fftw_complex***)q;
	
	fftw_plan* n_N = s->n_N;
	fftw_plan* m_M = s->m_M;
	fftw_plan* M_m = s->M_m;
	fftw_plan* E_e = s->E_e;
	fftw_plan** Q_q = s->Q_q;
	
	// store latest state in n and FFT
	for(u = 0; u < U; u++) {
		memcpy(n[u], m[u], Nxp*lNy*sizeof(double));
		fftw_execute(m_M[u]);
	}
	descale(s, m);
	
	// compute nonlocal components in Fourier space
	a = - 0.5/sigma/sigma;
	for(i = 0; i < lNx; i++) {
		gi = li0 + i;
		kx2 = gi*dkx;
		kx2 *= kx2;
		ij = Ny*i;
		for(j = 0; j < Ny; j++) {
			if(j < Ny/2) ky = j*dky;
			else ky = (j - Ny)*dky;
			k2 = kx2 + ky*ky;
			G = exp(a*k2);
			for(u = 0; u < U; u++) {
				// smoothed density
				E[u][ij] = G*M[u][ij];
				// gradient term coupling
				for(v = u + 1; v < U; v++) {
					d = mu2[u][v] - k2;
					Q[u][v][ij] = d*d*M[v][ij];
				}
				// multimode gradient term
				P = 1.0;
				for(w = 0; w < W[u]; w++) {
					p = nu2[u][w] - k2;
					P *= omega[u][w]*p*p + phi[u][w];
				}
				M[u][ij] *= P;
			}
			ij++;
		}
	}
	
	// FFT^(-1)
	for(u = 0; u < U; u++) {
		fftw_execute(M_m[u]);
		fftw_execute(E_e[u]);
		for(v = u + 1; v < U; v++) {
			fftw_execute(Q_q[u][v]);
		}
	}
	
	// zero average densities
	*fave = 0.0;
	memset(nave, 0.0, U*sizeof(double));
	
	// compute final result in direct space
	for(j = 0; j < lNy; j++) {
		ij = Nxp*j;
		for(i = 0; i < Nx; i++) {
			f[ij] = 0.0;
			for(u = 0; u < U; u++) {
				nunu = n[u][ij]*n[u][ij];
				f[ij] += lambda[u]*n[u][ij]
						+ 0.5*alpha[u][u]*nunu
						+ 0.5*beta[u][u]*n[u][ij]*m[u][ij]
						+ _3*gamma[u][u]*n[u][ij]*nunu
						+ 0.25*delta[u]*nunu*nunu;
				for(v = u + 1; v < U; v++) {
					nvnv = n[v][ij]*n[v][ij];
					f[ij] += alpha[u][v]*n[u][ij]*n[v][ij]
							+ beta[u][v]*n[u][ij]*q[u][v][ij]
							+ 0.5*gamma[u][v]*(nunu*n[v][ij] + n[u][ij]*nvnv)
							+ epsilon[u][v]*e[u][ij]*e[v][ij];
				}
				m[u][ij] = n[u][ij];
				nave[u] += n[u][ij];
			}
			*fave += f[ij];
			ij++;
		}
	}
	
	// communication between processes
	MPI_Allreduce(MPI_IN_PLACE, fave, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, nave, U, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	*fave *= _NxNy;
	for(u = 0; u < U; u++) {
		nave[u] *= _NxNy;
		fftw_execute(n_N[u]);	// restore solution in Fourier space
	}
	descale(s, n);
	
}

// updates solution in Fourier space according to the spectral method in Sec. B.4 in the book Phase-field methods in material science and engineering by Provatas and Elder (http://www.physics.mcgill.ca/~provatas/papers/Phase_Field_Methods_text.pdf)
void step(struct S* s) {

	int i, j, ij, u, v, gi;
	
	int t = s->t;
	int U = s->U;
	int Nx = s->Nx;
	int Ny = s->Ny;
	int Nxp = 2*(Nx/2 + 1);
	int lNx = s->lNx;
	int lNy = s->lNy;
	int li0 = s->li0;
	
	double m2, a, kx2, ky, k2, G2, d2;
	double NxNy = (double)Nx*Ny;
	double dx = s->dx;
	double dy = s->dy;
	double dkx = 2.0*pi/(dx*Nx);
	double dky = 2.0*pi/(dy*Ny);
	
	double sigma = s->sigma;
	double* lambda = s->lambda;
	double* delta = s->delta;
	double** alpha = s->alpha;
	double** beta = s->beta;
	double** gamma = s->gamma;
	double** epsilon = s->epsilon;
	double** mu2 = s->mu2;
	
	double* mc = s->mc;
	double** n = s->n;
	double** m = s->m;
	double** A = s->A;
	double** B = s->B;
	fftw_complex** N = (fftw_complex**)n;
	fftw_complex** M = (fftw_complex**)m;
	
	fftw_plan* n_N = s->n_N;
	fftw_plan* m_M = s->m_M;
	fftw_plan* M_m = s->M_m;

	// nan check
	if(m[0][0] != m[0][0]) fail(s, 39);
	
	// numerical instability fix
	if(t%100 == 0) {
		for(u = 0; u < U; u++) {
			memcpy(n[u], m[u], Nxp*lNy*sizeof(double));
			fftw_execute(n_N[u]);
		}
		descale(s, n);
	}
	
	// nonlinear part computed in direct space
	for(j = 0; j < lNy; j++) {
		ij = Nxp*j;
		for(i = 0; i < Nx; i++) {
			// reset array for couplings
			memset(mc, 0.0, U*sizeof(double));
			// compute couplings
			for(u = 0; u < U - 1; u++) {
				for(v = u + 1; v < U; v++) {
					mc[u] += alpha[u][v]*m[v][ij]
							+ 0.5*gamma[u][v]*(2.0*m[u][ij]*m[v][ij] + m[v][ij]*m[v][ij]);
					mc[v] += alpha[u][v]*m[u][ij]
							+ 0.5*gamma[u][v]*(2.0*m[v][ij]*m[u][ij] + m[u][ij]*m[u][ij]);
				}
			}
			// compute ideal part
			for(u = 0; u < U; u++) {
				m2 = m[u][ij]*m[u][ij];
				m[u][ij] = lambda[u] + gamma[u][u]*m2 + delta[u]*m[u][ij]*m2 + mc[u];
			}
			ij++;
		}
	}
	
	// FFT
	for(u = 0; u < U; u++) fftw_execute(m_M[u]);
	
	// update in Fourier space
	a = - 0.5/sigma/sigma;
	for(i = 0; i < lNx; i++) {
		gi = li0 + i;
		kx2 = gi*dkx;
		kx2 *= kx2;
		ij = Ny*i;
		for(j = 0; j < Ny; j++) {
			if(j < Ny/2) ky = j*dky;
			else ky = (j - Ny)*dky;
			k2 = kx2 + ky*ky;
			G2 = exp(a*k2);
			G2 *= G2;
			// compute couplings
			for(u = 0; u < U - 1; u++) {
				for(v = u + 1; v < U; v++) {
					d2 = mu2[u][v] - k2;
					d2 *= d2;
					M[u][ij] += (beta[u][v]*d2*N[v][ij] + epsilon[u][v]*G2*N[v][ij])*NxNy;
					M[v][ij] += (beta[u][v]*d2*N[u][ij] + epsilon[u][v]*G2*N[u][ij])*NxNy;
				}
			}
			// final result
			for(u = 0; u < U; u++) {
				N[u][ij] = A[u][ij]*N[u][ij] + B[u][ij]*M[u][ij];
				M[u][ij] = N[u][ij];	// copy result to M
			}
			ij++;
		}
	}
	
	// FFT^(-1)
	for(u = 0; u < U; u++) fftw_execute(M_m[u]);
	
}

// optimizes the system dimensions to eliminate strain
void optimize(struct S* s) {
	int xopt = s->xopt;	// optimized in x direction?
	int yopt = s->yopt;	// optimized in y direction?
	// sample sizes
	double dxs[] = {s->dx, s->dx - s->d, s->dx + s->d, s->dx, s->dx};
	double dys[] = {s->dy, s->dy, s->dy, s->dy - s->d, s->dy + s->d};
	double fs[5];							// average free energy densities of sizes sampled
	int i;
	for(i = 0; i < 5; i++) {				// sample sizes
		if(xopt == 0 && (i == 1 || i == 2)) continue;
		if(yopt == 0 && (i == 3 || i == 4)) continue;
		s->dx = dxs[i];
		s->dy = dys[i];
		evaluate(s);						// compute energy
		fs[i] = s->fave;					// save average free energy density
	}
	// interpolate new dx and dy (minimum of f = a+b*dx+c*dy+d*dx^2+e*dy^2)
	if(xopt != 0) s->dx = (s->d*(fs[1] - fs[2]) + 2.0*dxs[0]*(- 2.0*fs[0] + fs[1] +fs[2]))/(2.0*(fs[1] - 2.0*fs[0] + fs[2]));
	if(yopt != 0) s->dy = (s->d*(fs[3] - fs[4]) + 2.0*dys[0]*(- 2.0*fs[0] + fs[3] +fs[4]))/(2.0*(fs[3] - 2.0*fs[0] + fs[4]));
	// check that change in discretization is acceptable
	double l0 = 4.0*pi/sqrt(3.0);			// approximate dimensionless lenght scale (assuming hexagonal/honeycomb lattice and nu = 1)
	double dw = s->Nx*(s->dx - dxs[0]);		// change in horizontal system size
	double dh = s->Ny*(s->dy - dys[0]);		// ... vertical ...
	double dr = sqrt(dw*dw + dh*dh);			// "change vector"
	double x = 0.25*l0;						// limit to 1/4 of lattice constant (to ensure stability)
	if(dr > x) {							// if the change in system dimensions exceeds 1/4 of the lattice constant ...
		x /= dr;
		s->dx = x*s->dx + (1.0 - x)*dxs[0];	// ... truncate the change to 1/4 of the lattice constant
		s->dy = x*s->dy + (1.0 - x)*dys[0];
	}
	// update A and B
	update_AB(s);							// dx and dy changed -> need to update
	// update sampling step size (tries to keep it in the same ballpark with dr (for (hopefully) more accurate optimization))
	double ddx = s->dx - dxs[0];
	double ddy = s->dy - dys[0];
	double ddr = sqrt(ddx*ddx + ddy*ddy);	// discretization change
	if(ddr < s->d) s->d *= 0.5;				// if change vector < d, halve d
	else s->d *= 2.0;						// otherwise double
	if(s->d < 1.0e-6) s->d *= 2.0;			// can cause numerical issues if d gets too small
}

// relaxes the system
void relax(struct S* s) {
	int n = 0;
	// initial optimization sampling step
	s->d = 0.0001;
	for(s->t = 0; s->t <= s->tend; s->t++) {
		// optimize?
		if(s->Topt > 0 && s->t > 0 && s->t%s->Topt == 0) optimize(s);
		// print output?
		if(s->t%s->Tprint == 0) {
			evaluate(s);
			print_output(s);
		}
		// write out data?
		if(s->t == s->twrite[n]) {
			write_data(s);
			if(s->fname != NULL) write_energy(s);
			n++;
		}
		if(s->t < s->tend) step(s);
	}
}

// frees memory
void clean(struct S* s) {

	int u, v;
	int U = s->U;
	
	fftw_free(s->twrite);
	
	fftw_free(s->mc);
	fftw_free(s->nave);
	
	fftw_free(s->W);
	fftw_free(s->delta);
	fftw_free(s->lambda);
	fftw_free(s->alpha[0]);
	fftw_free(s->beta[0]);
	fftw_free(s->gamma[0]);
	fftw_free(s->epsilon[0]);
	fftw_free(s->mu2[0]);
	fftw_free(s->alpha);
	fftw_free(s->beta);
	fftw_free(s->gamma);
	fftw_free(s->epsilon);
	fftw_free(s->mu2);
	
	for(u = 0; u < U; u++) {
		fftw_free(s->innames[u]);
		fftw_free(s->outnames[u]);
	
		fftw_free(s->nu2[u]);
		fftw_free(s->omega[u]);
		fftw_free(s->phi[u]);
		
		fftw_free(s->A[u]);
		fftw_free(s->B[u]);
		fftw_free(s->n[u]);
		fftw_free(s->m[u]);
		fftw_free(s->e[u]);
		
		for(v = u + 1; v < U; v++) {
			fftw_free(s->q[u][v]);
			fftw_destroy_plan(s->Q_q[u][v]);
		}
		
		fftw_free(s->q[u]);
	
		fftw_destroy_plan(s->n_N[u]);
		fftw_destroy_plan(s->m_M[u]);
		fftw_destroy_plan(s->M_m[u]);
		fftw_destroy_plan(s->E_e[u]);
		fftw_free(s->Q_q[u]);
	}
	
	fftw_free(s->fname);
	fftw_free(s->innames);
	fftw_free(s->outnames);
	
	fftw_free(s->nu2);
	fftw_free(s->omega);
	fftw_free(s->phi);
	
	fftw_free(s->f);
	fftw_free(s->A);
	fftw_free(s->B);
	fftw_free(s->n);
	fftw_free(s->m);
	fftw_free(s->e);
	fftw_free(s->q);

	fftw_free(s->n_N);
	fftw_free(s->m_M);
	fftw_free(s->M_m);
	fftw_free(s->E_e);
	fftw_free(s->Q_q);
	
}

int main(int argc, char** argv) {

	// initialize MPI
	MPI_Init(&argc, &argv);
	fftw_mpi_init();
	
	// save processes' ranks
	struct S s;
	MPI_Comm_rank(MPI_COMM_WORLD, &s.r);

	// print instructions
	if(argc == 1) {
		if(s.r == 0) printf("This is a program for relaxing PFC systems. It expects a single argument specifying an input file. An input file in turn specifies the input and output of a relaxation, model parameters and relaxation parameters. This program supports multiple length scales and multicomponent systems with multiple density fields and different couplings between them, including a new coupling that enables straightforward modeling of heterostructures (see arXiv:1908.05564). Both conserved and nonconserved dynamics are supported. See the source code for definitions of the model and the dynamics. The program also implements a system size optimization algorithm; see Sec. 4.3.1 in my thesis (available: http://urn.fi/URN:ISBN:978-952-60-8608-8). This program uses MPI and FFTW's fast Fourier transform routines for efficient computation.\n");
		MPI_Finalize();
		return 0;
	}
	
	// invalid syntax
	else if(argc != 2) {
		if(s.r == 0) printf("Error: Invalid syntax.\n");
		MPI_Abort(MPI_COMM_WORLD, 41);
	}
	
	// initialize some stuff
	s.t0 = time(NULL);						// save wall-clock time
	MPI_Comm_size(MPI_COMM_WORLD, &s.R);	// save number of processes
	strcpy(s.input, argv[1]);
	
	// read input file
	read_input(&s);
	
	// initialize
	initialize(&s);
	
	// relax
	relax(&s);
	
	// clean up
	clean(&s);
	MPI_Finalize();
	
	return 0;
	
}
