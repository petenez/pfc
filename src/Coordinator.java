// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package atomistic;

import java.io.*;
import java.util.*;

import static java.lang.System.exit;

//implements a "coordinator" for converting a PFC density field into atom coordinates
public class Coordinator {

	private static int[] u = new int[] {1, 1, 0, -1, -1, -1, 0, 1};	// arrays for looping through ...
	private static int[] v = new int[] {0, 1, 1, 1, 0, -1, -1, -1};	// ... neighboring indices

	private int Nx;		// PFC grid dimensions
	private int Ny;
	private int Bx;		// numbers of local minima bins
	private int By;
	private double w;	// system dimensions
	private double h;
	private double lA;	// lattice constant
	private double ux;	// length units
	private double uy;
	private ArrayList<Atom> atoms = new ArrayList<Atom>();	// list of atoms
	private ArrayList<Atom> rings = new ArrayList<Atom>();	// list of rings
	private ArrayList<ArrayList<ArrayList<Atom>>> bins;		// local minima bins

	public Coordinator(int Nx, int Ny, double dx, double dy, double lPFC, double lA) {
		this.Nx = Nx;
		this.Ny = Ny;
		this.lA = lA;
		ux = dx*lA/lPFC;	// grid point width in the units of lA
		uy = dy*lA/lPFC;
		w = Nx*ux;
		h = Ny*uy;
		this.Bx = (int)Math.ceil(0.5*w/lA);	// bin size ~2*lA
		this.By = (int)Math.ceil(0.5*h/lA);
		if(Bx == 1 || By == 1) {
			System.out.println("Warning: The conversion may fail for such small system sizes.");
		}

		bins = new ArrayList<ArrayList<ArrayList<Atom>>>(By);	// initialize bins
		for(int m = 0; m < By; m++) {
			bins.add(new ArrayList<ArrayList<Atom>>(Bx));
			for(int n = 0; n < Bx; n++) {
				bins.get(m).add(new ArrayList<Atom>());
			}
		}
	}

	// wraps x-coordinates
	public int wx(int x) {
		if(x < 0) return x + Nx;
		if(x >= Nx) return x - Nx;
		return x;
	}

	// wraps y-coordinates
	public int wy(int y) {
		if(y < 0) return y + Ny;
		if(y >= Ny) return y - Ny;
		return y;
	}

	// wraps bin x-indices
	public int wnx(int nx) {
		if(nx < 0) return nx + Bx;
		if(nx >= Bx) return nx - Bx;
		return nx;
	}

	// wraps bin y-indices
	public int wny(int ny) {
		if(ny < 0) return ny + By;
		if(ny >= By) return ny - By;
		return ny;
	}

	// wraps a vector wrt. system dimensions
	public double[] w(double[] r) {
		double x = r[0];
		double y = r[1];
		if(x < 0.0) x += w;
		else if(x >= w) x -= w;
		if(y < 0.0) y += h;
		else if(y >= h) y -= h;
		return new double[] {x, y};
	}

	// returns the bin x-index corresponding to vector r
	public int nx(double[] r) {
		return (int)(r[0]/w*Bx);
	}

	// returns the bin y-index corresponding to vector r
	public int ny(double[] r) {
		return (int)(r[1]/h*By);
	}

	// returns the bin indexed by nx and ny
	public ArrayList<Atom> bin(int nx, int ny) {
		return bins.get(ny).get(nx);
	}

	// returns the bin corresponding to vector r
	public ArrayList<Atom> bin(double[] r) {
		return bins.get((int)ny(r)).get((int)nx(r));
	}

	// computes the difference vector from a to b taking periodic boundary conditions into account
	public double[] diff(double[] b, double[] a) {
		double[] ab = V.sub(b, a);
		if(ab[0] > 0.5*w) ab[0] -= w;
		else if(ab[0] < -0.5*w) ab[0] += w;
		if(ab[1] > 0.5*h) ab[1] -= h;
		else if(ab[1] < -0.5*h) ab[1] += h;
		return ab;
	}

	// adds nth as a member of the three rings, adds the three rings as neighbors if not previously neighbors
	public void bros_n_hos(Atom nth, Atom ith, Atom jth, Atom kth) {
		nth.add_ho(ith);
		nth.add_ho(jth);
		nth.add_ho(kth);
		ith.add_ho(nth);
		jth.add_ho(nth);
		kth.add_ho(nth);
		if(!ith.is_bro(jth)) {
			ith.add_bro(jth);
			jth.add_bro(ith);
		}
		if(!ith.is_bro(kth)) {
			ith.add_bro(kth);
			kth.add_bro(ith);
		}
		if(!jth.is_bro(kth)) {
			jth.add_bro(kth);
			kth.add_bro(jth);
		}
	}

	// Loads the PFC system, detects local minima, finds nearest-neighbor triplets of them, places atoms in the center of each triplet, finds their neighbors and checks that there's a reasonable number of neighbors.
	public void load(String input) throws IOException {
		double[][] psi = new double[Ny][Nx];
		BufferedReader reader = new BufferedReader(new FileReader(input));
		reader.readLine();
		for(int j = 0; j < Ny; j++) {		// load PFC density data
			for(int i = 0; i < Nx; i++) {
				psi[j][i] = Double.parseDouble(reader.readLine());
			}
		}
		reader.close();
		Atom nth;
		for(int j = 0; j < Ny; j++) {
			outer:
			for(int i = 0; i < Nx; i++) {
				for(int n = 0; n < 8; n++) {
					// if a neighbor's density is lower, not a local minimum
					if(psi[wy(j + v[n])][wx(i + u[n])] < psi[j][i]) {
						continue outer;
					}
				}

				// use quadratic interpolation to estimate more accurate location
				double x0 = i - 1;	// neighbor
				double x1 = i;	// current minimum
				double x2 = i + 1;	// neighbor
				double x02 = x0*x0;
				double x12 = x1*x1;
				double x22 = x2*x2;
				double p0 = psi[j][wx(i - 1)];	// periodic boundary conditions
				double p1 = psi[j][i];
				double p2 = psi[j][wx(i + 1)];
				double A = (x2*(p1 - p0) + x1*(p0 - p2) + x0*(p2 - p1));	// y = A*x^2 + B*x + C
				double B = (x22*(p0 - p1) + x12*(p2 - p0) + x02*(p1 - p2));
				double x = - B/(2.0*A)*ux;	// extremum at x = -B/(2A)
				if(Math.abs(A) < 1.0e-14) x = i;	// linear fit

				x0 = j - 1;	// same in the other direction
				x1 = j;
				x2 = j + 1;
				x02 = x0*x0;
				x12 = x1*x1;
				x22 = x2*x2;
				p0 = psi[wy(j - 1)][i];
				p1 = psi[j][i];
				p2 = psi[wy(j + 1)][i];
				A = (x2*(p1 - p0) + x1*(p0 - p2) + x0*(p2 - p1));
				B = (x22*(p0 - p1) + x12*(p2 - p0) + x02*(p1 - p2));
				double y = - B/(2.0*A)*uy;
				if(Math.abs(A) < 1.0e-14) y = j;
				nth = new Atom(rings.size(), w(new double[] {x, y}));
				rings.add(nth);		// add local minimum to corresponding list
				bin(nth.r()).add(nth);	// add to corresponding bin

				// density data might have limited precision whereby neighboring grid points might have identical values, set current density and those of neighbors very low to ensure uniqueness of current local minimum
				psi[j][i] -= 2000.0;
				for(int n = 0; n < 8; n++) {
					psi[wy(j + v[n])][wx(i + u[n])] -= 1000.0;
				}
			}
		}

		// Find triplets of local minima that are closest neighbors to each other. Practically find the triple points of three local minima (where each one is equally distant) and make sure that no fourth atom atom is closer to it. Place a atom atom at the average of the three minimas' coordinates. Find the neighbors of the atom atoms and minima (which actually correspond to atom atom rings) and check that the bonding makes sense.
		Atom ith, jth, kth, lth;	// 'Points' for minima and atoms
		int nx, ny;								// indices for bins
		double D, _D, Ux, Uy, R2;
		double[] ij, ik, it, r;			// vectors
		ArrayList<Atom> binj, bink, binl;	// bins
		for(int i = 0; i < rings.size(); i++) {	// go through rings
			ith = rings.get(i);
			nx = nx(ith.r());
			ny = ny(ith.r());
			for(int vj = -1; vj < 2; vj++) {	// go through its surrounding bins
				for(int uj = -1; uj < 2; uj++) {
					binj = bin(wnx(nx + uj), wny(ny + vj));
					for(int j = 0; j < binj.size(); j++) {	// go through minima in current bin
						jth = binj.get(j);
						if(jth.k() <= ith.k()) continue;	// consider each triplet only once
						// compute some vectors to find the triple point
						ij = diff(jth.r(), ith.r());
						if(V.abs(ij) > 2.0*lA) continue;	// closest neighbors can't be farther apart
						for(int vk = -1; vk < 2; vk++) {	// again check bins
							for(int uk = -1; uk < 2; uk++) {
								bink = bin(wnx(nx + uk), wny(ny + vk));
								k_loop:
								// rings in this bin
								for(int k = 0; k < bink.size(); k++) {
									kth = bink.get(k);
									if(kth.k() <= jth.k()) continue;
									// compute vectors
									ik = diff(kth.r(), ith.r());
									if(V.abs(ik) > 2.0*lA) continue;
									if(V.abs(diff(kth.r(), jth.r())) > 2.0*lA) continue;
									// for expressions below see: https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_2
									D = 2.0*(ij[0]*ik[1] - ij[1]*ik[0]);
									if(D == 0.0) continue;	// the three minima are lined up
									_D = 1.0/D;
									Ux = (ik[1]*(ij[0]*ij[0] + ij[1]*ij[1]) - ij[1]*(ik[0]*ik[0] + ik[1]*ik[1]))*_D;
									Uy = (ij[0]*(ik[0]*ik[0] + ik[1]*ik[1]) - ik[0]*(ij[0]*ij[0] + ij[1]*ij[1]))*_D;
									it = new double[] {Ux, Uy};	// position of triple point (from ith minimum)

									// check that no fourth minimum is closer to the triple point
									R2 = V.abs2(it);
									for(int vl = -1; vl < 2; vl++) {
										for(int ul = -1; ul < 2; ul++) {
											binl = bin(wnx(nx + ul), wny(ny + vl));
											for(int l = 0; l < binl.size(); l++) {
												lth = binl.get(l);
												if(lth == ith || lth == jth || lth == kth) continue;
												if(V.abs2(V.sub(diff(lth.r(), ith.r()), it)) < R2) continue k_loop;
											}
										}
									}

									// add atom atom, update neighbors
									nth = new Atom(atoms.size(), w(V.sum(ith.r(), V.mul(1.0/3.0, V.sum(ij, ik)))));
									atoms.add(nth);
									bros_n_hos(nth, ith, jth, kth);
								}
							}
						}
					}
				}
			}
		}

		// find neighboring atoms of each atom
		Atom mth;
		for(int n = 0; n < atoms.size(); n++) {	// all atoms
			nth = atoms.get(n);
			for(int i = 0; i < nth.count_hos(); i++) {	// ith ring its a member of
				ith = nth.get_ho(i);
				for(int j = i + 1; j < nth.count_hos(); j++) {	// jth ring its a member of
					jth = nth.get_ho(j);
					for(int k = 0; k < ith.count_hos(); k++) {	// kth atom of ith ring
						mth = ith.get_ho(k);
						// if its not C, but shares two rings with it and is not yet a neighbor ...
						if(nth != mth && mth.is_ho(jth) && !nth.is_bro(mth)) {
							nth.add_bro(mth);	// ... add as neighbors
							mth.add_bro(nth);
						}
					}
				}
			}
		}

		// check if every atom has three neighbors and if each ring has 5-7 neighboring rings etc.
		boolean flag1 = false;
		boolean flag2 = false;
		boolean flag3 = false;
		for(int n = 0; n < atoms.size(); n++) {
			nth = atoms.get(n);
			if(!flag1 && nth.count_bros() != 3) {
				System.out.println("Warning: An atom with fewer or more than three bonds was detected.");
				flag1 = true;
			}
			if(!flag2 && nth.count_hos() != 3) {
				System.out.println("Warning: An atom shared by fewer or more than three rings was detected.");
				flag2 = true;
			}
		}
		for(int i = 0; i < rings.size(); i++) {
			ith = rings.get(i);
			if(!flag3 && (ith.count_bros() < 5 || ith.count_bros() > 7)) {
				System.out.println("Warning: A ring with fewer than five or more than seven atoms was detected.");
				flag3 = true;
			}
			if(ith.count_bros() != ith.count_hos()) {
				System.out.println("Error: Mismatching numbers of atoms and neighboring rings were detected for ring " + ith.k() + ".");
				exit(1);
			}
		}
	}

	// writes out the atom coordinates
	// format: atom type, its coordinates and indices of its neighbors
	public void write_atoms(String output) throws IOException {
		Atom nth, mth;
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		writer.write(w + " " + h + "\n");
		for(int n = 0; n < atoms.size(); n++) {
			nth = atoms.get(n);
			writer.write("A " + nth.r()[0] + " " + nth.r()[1] + " 0.0");
			for(int m = 0; m < nth.count_bros(); m++) {
				mth = nth.get_bro(m);
				writer.write(" " + mth.k());
			}
			writer.write("\n");
		}
		writer.close();
	}

	// writes out the nonhexagon coordinates
	// format: indices of the member atoms
	public void write_nonhexagons(String output) throws IOException {
		Atom ith, nth;
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		for(int i = 0; i < rings.size(); i++) {
			ith = rings.get(i);
			if(ith.count_bros() == 6) continue;
			for(int n = 0; n < ith.count_hos(); n++) {
				nth = ith.get_ho(n);
				writer.write(nth.k() + " ");
			}
			writer.write("\n");
		}
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if(args.length == 0) {
			System.out.println("This is a tool for converting PFC density fields into atom coordinates. This tool can be applied when honeycomb structures are modeled using inverted hexagonal density fields (their maxima form a honeycomb lattice); see Sec. 4.3.1 in my thesis (available: http://urn.fi/URN:ISBN:978-952-60-8608-8).");
			System.out.println();
			System.out.println("Examples of valid syntax:");
			System.out.println();
			System.out.println("java -jar coordinator.jar density.n coordinates.xyz [input length scale] [output length scale]");
			System.out.println();
			System.out.println("java -jar coordinator.jar density.n coordinates.xyz nonhexagons.nh [input length scale] [output length scale]");
			System.out.println();
			System.out.println("The program expects the names of the input density and the output coordinate files as its first arguments. A coordinate file lists each atom's type (the program assigns all atoms with a label \"A\"), its coordinates and indices of its neighbors. The program can also write out the nonhexagonal atom rings (using the indices of the member atoms) if a filename is specified (the second example above). The input length scale is the length scale of the PFC structures, typically 4*pi/sqrt(3) ~ 7.3. The output length scale is the lattice constant in the units desired, e.g., approx. 2.5 Å or 0.25 nm for graphene.");
			exit(0);
		}
		if(args.length != 4 && args.length != 5) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

		int k = 0;
		String input = args[k++];	// input file
		File file = new File(input);
		if(!file.exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}

		BufferedReader reader = new BufferedReader(new FileReader(input));
		// reads dimensions and discretization from first line
		String[] temp = reader.readLine().trim().split("\\s+");
		int Nx = Integer.parseInt(temp[0]);	// grid points in x direction
		int Ny = Integer.parseInt(temp[1]);	// ... y ...
		double dx = Double.parseDouble(temp[2]);	// x-discretization
		double dy = Double.parseDouble(temp[3]);	// y...
		reader.close();

		String cout = args[k++];	// output file for atoms
		String nhout = null;
		if(args.length == 5) nhout = args[k++];	// output file for nonhexagonal rings

		double lPFC = Double.parseDouble(args[k++]);	// PFC length scale
		double lA = Double.parseDouble(args[k++]);	// lattice constant

		Coordinator coordinator = new Coordinator(Nx, Ny, dx, dy, lPFC, lA);
		coordinator.load(input);	// loads input data and processes it
		coordinator.write_atoms(cout);		// writes out atoms
		if(args.length == 5) coordinator.write_nonhexagons(nhout);	// writes out nonhexagonal rings
	}
}