// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package data;

import java.io.*;
import java.util.Random;

import static java.lang.System.exit;

public class PFCInit {

	public static int s2i(String s) {
		int i = 0;
		try {
			i = Integer.parseInt(s);
		} catch (NumberFormatException e) {
			System.out.println("Error: Invalid input file. Integer expected instead of \"" + s + "\".");
			exit(1);
		}
		return i;
	}

	public static double s2d(String s) {
		double d = 0.0;
		try {
			d = Double.parseDouble(s);
		} catch (NumberFormatException e) {
			System.out.println("Error: Invalid input file. Decimal number expected instead of \"" + s + "\".");
			exit(1);
		}
		return d;
	}

	public static int w(int i, int N) {
		while(i < 0) i += N;
		while(i >= N) i -= N;
		return i;
	}

	public static double w(double d, double L) {
		while(d < 0) d += L;
		while(d >= L) d -= L;
		return d;
	}

	public static double OMA(double x, double y, double sym, double a) {
		double q = 2.0*Math.PI/a;
		if(sym == 3) {
			double qx = q;
			double qy = q/Math.sqrt(3.0);
			return Math.cos(qx*x)*Math.cos(qy*y) + 0.5*Math.cos(2.0*qy*y);
		}
		else if(sym == 4) {
			return Math.cos(q*x)*Math.cos(q*y);
		}
		else if(sym == 6) {
			double qx = q;
			double qy = q/Math.sqrt(3.0);
			return - 2.0*Math.cos(qx*x - qy*y) - 2.0*Math.cos(2.0*qy*y) - 2.0*Math.cos(- qx*x - qy*y) + Math.cos(- qx*x - 3.0*qy*y) + Math.cos(2.0*qx*x) + Math.cos(- qx*x + 3.0*qy*y);
		}
		else {
			System.out.println("Error: Invalid lattice symmetry specified. Only triangular, square and hexagonal symmetries are supported.");
			exit(1);
		}
		return 0.0;
	}

	public static int closest(double u, double v, double[] x, double[] y, double W, double H) {
		int N = x.length;
		int n = 0;
		double dx = u - x[0];
		if(dx > 0.5*W) dx -= W;
		else if(dx < -0.5*W) dx += W;
		double dy = v - y[0];
		if(dy > 0.5*H) dy -= H;
		else if(dy < -0.5*H) dy += H;
		double d2 = dx*dx + dy*dy;
		double d2_;
		for(int m = 1; m < N; m++) {
			dx = u - x[m];
			if(dx > 0.5*W) dx -= W;
			else if(dx < -0.5*W) dx += W;
			dy = v - y[m];
			if(dy > 0.5*H) dy -= H;
			else if(dy < -0.5*H) dy += H;
			d2_ = dx*dx + dy*dy;
			if(d2_ < d2) {
				n = m;
				d2 = d2_;
			}
		}
		return n;
	}

	public static void main(String[] args) throws IOException {

		if(args.length != 1) {
			System.out.println("This is a tool for constructing 2D density fields for PFC calculations. The program expects one argument specifying an input file. The input file is expected to begin with a specification of system dimensions and spatial discretization, and to end with a specification of the output file and its decimal precision. Between these, multiple grains can be constructed and the seed number for the pseudo random number generator can be specified. Empty and commented (#) lines are allowed.");
			System.out.println();
			System.out.println("An example of a simple input file:");
			System.out.println();
			System.out.println(
				"# system dimensions and discretization\n" +
				"D\n" +
				"\n" +
				"#\tNx\t\tNy\t\tdx\t\tdy\n" +
				"\t512\t\t512\t\t0.7\t\t0.7\n" +
				"\n" +
				"# RNG seed\n" +
				"S\n" +
				"\n" +
				"#\tseed number\n" +
				"\t12345\n" +
				"\n" +
				"# circular grain\n" +
				"C\n" +
				"\n" +
				"#\tx\t\ty\t\tR\t\ttype\n" +
				"\t0.0\t\t0.0\t\t100.0\t0\n" +
				"\t\n" +
				"#\tsym\t\ta\t\tnave\tamp\t\tnamp\t\n" +
				"\t3\t\t7.2552\t0.0\t\t1.0\t\t0.1\n" +
				"\t\n" +
				"#\ttx1\t\tty1\t\ttheta\txr\t\tyr\t\ttx2\t\tty2\n" +
				"\t0.0\t\t0.0\t\t0.0\t\t0.0\t\t0.0\t\t0.0\t\t0.0\n" +
				"\n" +
				"# output\n" +
				"O\n" +
				"\n" +
				"#\tfilename\tprecision\n" +
				"\texample.n\t6");
			System.out.println();
			System.out.println("Label D denotes specification of the system dimensions and spatial discretization. Integers Nx and Ny give the system dimensions in grid points, and decimal numbers dx and dy give the spatial discretization in x and y directions.");
			System.out.println();
			System.out.println("Label S denotes specification of the seed number for the pseudo random number generator (java.util.Random). If left unspecified, the generator is initialized without a seed number. The generator is needed if white noise is to be added to the grains.");
			System.out.println();
			System.out.println("Label C denotes specification of a circular grain. Decimal numbers x, y and R give the x and y coordinates of the center point and the radius (note that dimensionless length units are expected). Integers type and sym give the write type and the lattice symmetry of the grain. Type can be 0, 1 or 2, corresponding to overwriting any previous structures coinciding with the new grain, to writing the new grain in superposition with any coinciding previous structures and to writing the new grain as a phantom one leaving any coinciding previous structures unchanged (can be useful with Voronoi initialization; see below). Supported lattice symmetries are triangular or hexagonal (3), square (4) and honeycomb (6). Decimal numbers a, nave, amp and namp give the length scale, average density, amplitude and noise amplitude of the grain. Decimal numbers tx1, ty1, theta, xr, yr, tx2 and ty2 give the x and y translation of the lattice structure before rotation specified by the angle (degrees) and x and y coordinates of the rotation axis, and the x and y translation after the rotation.");
			System.out.println();
			System.out.println("Label O denotes specification of the output file. Filename and decimal precision are expected.");
			System.out.println();
			System.out.println("A second example:");
			System.out.println();
			System.out.println(
				"# system dimensions and discretization\n" +
				"D\n" +
				"\n" +
				"#\tNx\t\tNy\t\tdx\t\tdy\n" +
				"\t512\t\t512\t\t0.7\t\t0.7\n" +
				"\n" +
				"# load saved system\n" +
				"L\n" +
				"\n" +
				"#\tfilename\ttype\n" +
				"\texample.n\t\t0\n" +
				"\t\n" +
				"#\ti0\t\tj0\t\tdNx\t\tdNy\t\ti1\t\tj1\n" +
				"\t-150\t-100\t200\t\t200\t\t100\t\t50\n" +
				"\n" +
				"# transform\n" +
				"T\n" +
				"\n" +
				"#\tnave\tfamp\n" +
				"\t-0.5\t0.5\n" +
				"\t\n" +
				"# polygon grain\n" +
				"P\t4\n" +
				"\n" +
				"#\tx\t\ty\n" +
				"\t50.0\t250.0\n" +
				"\t300.0\t290.0\n" +
				"\t160.0\t350.0\n" +
				"\t0.0\t\t370.0\n" +
				"\t\n" +
				"#\ttype\n" +
				"\t1\t\n" +
				"\t\n" +
				"#\tsym\t\ta\t\tnave\tamp\t\tnamp\t\n" +
				"\t6\t\t7.2552\t-1.0\t1.0\t\t0.0\n" +
				"\t\n" +
				"#\ttx1\t\tty1\t\ttheta\txr\t\tyr\t\ttx2\t\tty2\n" +
				"\t0.0\t\t0.0\t\t26.0\t0.0\t\t0.0\t\t0.0\t\t0.0\n" +
				"\n" +
				"# output\n" +
				"O\n" +
				"\n" +
				"#\tfilename\tprecision\n" +
				"\texample2.n\t6");
			System.out.println();
			System.out.println("This example demonstrates loading a part of the density field constructed in the previous example, placing it in the new density field and finally superposing with a polygon grain specified by its corners.");
			System.out.println();
			System.out.println("Label L denotes specification of the previous density field to be loaded. The part of the previous density field to be loaded is specified by the i and j indices of its bottom-left grid point (i0, j0), and by the part's dimensions dNx x dNy. The part is placed into the new density field so that its bottom-left grid point is (i1, j1).");
			System.out.println();
			System.out.println("Label T denotes specification of transforming the current density field. Decimal numbers nave and famp give the new average density and the scaling factor for the amplitude.");
			System.out.println();
			System.out.println("Label P denotes specification of a polygon grain given by its (in this case 4) corner points (x and y coordinates). Note that here type = 1, meaning that the new grain will be written in superposition with any previous structures.");
			System.out.println();
			System.out.println("A third example:");
			System.out.println();
			System.out.println(
				"# system dimensions and discretization\n" +
				"D\n" +
				"\n" +
				"#\tNx\t\tNy\t\tdx\t\tdy\n" +
				"\t512\t\t512\t\t0.7\t\t0.7\n" +
				"\n" +
				"# load saved system\n" +
				"L\n" +
				"\n" +
				"#\tfilename\ttype\n" +
				"\texample2.n\t\t0\n" +
				"\t\n" +
				"#\ti0\t\tj0\t\tdNx\t\tdNy\t\ti1\t\tj1\n" +
				"\t0\t\t0\t\t512\t\t512\t\t0\t\t0\n" +
				"\n" +
				"# Voronoi grains (3)\n" +
				"V\t3\n" +
				"\n" +
				"# 1st grain\n" +
				"#\tx\t\ty\t\ttype\n" +
				"\t100.0\t100.0\t0\n" +
				"\t\n" +
				"#\tsym\t\ta\t\tnave\tamp\t\tnamp\t\n" +
				"\t3\t\t7.2552\t0.0\t\t1.0\t\t0.2\n" +
				"\t\n" +
				"#\ttx1\t\tty1\t\ttheta\txr\t\tyr\t\ttx2\t\tty2\n" +
				"\t0.0\t\t0.0\t\t5.0\t\t0.0\t\t0.0\t\t0.0\t\t0.0\n" +
				"\n" +
				"# 2nd grain\n" +
				"#\tx\t\ty\t\ttype\n" +
				"\t250.0\t400.0\t1\n" +
				"\t\n" +
				"#\tsym\t\ta\t\tnave\tamp\t\tnamp\t\n" +
				"\t3\t\t7.2552\t0.1\t\t1.0\t\t0.0\n" +
				"\t\n" +
				"#\ttx1\t\tty1\t\ttheta\txr\t\tyr\t\ttx2\t\tty2\n" +
				"\t0.0\t\t0.0\t\t10.0\t0.0\t\t0.0\t\t0.0\t\t0.0\n" +
				"\n" +
				"# 3rd grain\n" +
				"#\tx\t\ty\t\ttype\n" +
				"\t50.0\t300.0\t2\n" +
				"\t\n" +
				"#\tsym\t\ta\t\tnave\tamp\t\tnamp\t\n" +
				"\t3\t\t7.2552\t0.3\t\t1.0\t\t0.0\n" +
				"\t\n" +
				"#\ttx1\t\tty1\t\ttheta\txr\t\tyr\t\ttx2\t\tty2\n" +
				"\t0.0\t\t0.0\t\t25.0\t0.0\t\t0.0\t\t0.0\t\t0.0\n" +
				"\n" +
				"# output\n" +
				"O\n" +
				"\n" +
				"#\tfilename\tprecision\n" +
				"\texample3.n\t6");
			System.out.println();
			System.out.println("This example demonstrates a Voronoi initialization (denoted by label V) where a Voronoi diagram is formed to partition the system into (in this case 3) grains based on the seed points given to the grains. Note that different write types are used to demonstrate overwriting, superposing and leaving the grains as phantoms.");
			exit(0);
		}

		String name = args[0];
		File file = new File(name);
		if(!file.exists()) {
			System.out.println("Error: Input file \"" + name + "\" not found.");
			exit(1);
		}
		String line;
		String[] words;
		String[] steps = new String[100];
		String[] steps2;
		BufferedReader reader = new BufferedReader(new FileReader(name));

		// clean input
		int k = 0;
		outer:
		while((line = reader.readLine()) != null) {
			words = line.trim().split("\\s+");
			for(int l = 0; l < words.length; l++) {
				// empty
				if(words[l].isEmpty()) continue;
				// comment
				if(words[l].charAt(0) == '#') continue outer;
				steps[k] = words[l];
				k++;
				if(k == steps.length) {
					steps2 = new String[10*steps.length];
					for(int m = 0; m < k; m++) {
						steps2[m] = steps[m];
					}
					steps = steps2;
				}
			}
		}
		int K = k;

		// build system
		Random rng = new Random();
		int Nx = -1;
		int Ny = -1;
		double dx = 0.0;
		double dy = 0.0;
		double[][] arr = null;
		k = 0;
		while(k < K) {

			// dimensions and discretization
			if(steps[k].equals("D")) {
				k++;
				Nx = s2i(steps[k]);
				k++;
				Ny = s2i(steps[k]);
				k++;
				dx = s2d(steps[k]);
				k++;
				dy = s2d(steps[k]);
				k++;
				arr = new double[Ny][Nx];
			}

			// seed
			else if(steps[k].equals("S")) {
				k++;
				rng = new Random(s2i(steps[k]));
				k++;
			}

			// load saved
			else if(steps[k].equals("L")) {
				k++;
				name = steps[k];
				k++;
				int type = s2i(steps[k]);
				k++;
				if(type == 2) {
					System.out.println("Warning: Phantom type used for the system from data file \"" + name + "\".");
					continue;
				}
				if(type < 0 || type > 2) {
					System.out.println("Warning: Invalid type used for the system from data file \"" + name + "\". Defaulting to overwrite.");
				}

				file = new File(name);
				if(!file.exists()) {
					System.out.println("Error: Data file \"" + name + "\" not found.");
					exit(1);
				}
				BufferedReader reader2 = new BufferedReader(new FileReader(name));
				line = reader2.readLine();
				words = line.split("\\s+");
				int Nx2 = s2i(words[0]);
				int Ny2 = s2i(words[1]);
				double dx2 = s2d(words[2]);
				double dy2 = s2d(words[3]);
				if(dx2 != dx || dy2 != dy) {
					System.out.println("Warning: System from data file \"" + name + "\" has incompatible discretization.");
				}

				int i0 = w(s2i(steps[k]), Nx2);
				k++;
				int j0 = w(s2i(steps[k]), Ny2);
				k++;
				int dNx = s2i(steps[k]);
				k++;
				int dNy = s2i(steps[k]);
				k++;
				int i1 = w(s2i(steps[k]), Nx);
				k++;
				int j1 = w(s2i(steps[k]), Ny);
				k++;
				int is, js;
				double[][] arr2 = new double[Ny2][Nx2];
				for(int j = 0; j < Ny2; j++) {
					for(int i = 0; i < Nx2; i++) {
						arr2[j][i] = s2d(reader2.readLine());
					}
				}
				int it, jt;
				for(int j = 0; j < dNy; j++) {
					js = w(j0 + j, Ny2);
					jt = w(j1 + j, Ny);
					for(int i = 0; i < dNx; i++) {
						is = w(i0 + i, Nx2);
						it = w(i1 + i, Nx);
						arr[jt][it] = type*arr[jt][it] + arr2[js][is];
					}
				}
				reader2.close();

			}

			// Voronoi grains
			else if(steps[k].equals("V")) {
				k++;
				int N = s2i(steps[k]);
				k++;
				double[] x = new double[N];
				double[] y = new double[N];
				int[] type = new int[N];
				int[] sym = new int[N];
				double[] a = new double[N];
				double[] nave = new double[N];
				double[] amp = new double[N];
				double[] namp = new double[N];
				double[] tx1 = new double[N];
				double[] ty1 = new double[N];
				double[] theta = new double[N];
				double[] xr = new double[N];
				double[] yr = new double[N];
				double[] tx2 = new double[N];
				double[] ty2 = new double[N];
				double W = Nx*dx;
				double H = Ny*dy;
				for(int n = 0; n < N; n++) {
					x[n] = w(s2d(steps[k]), W);
					k++;
					y[n] = w(s2d(steps[k]), H);
					k++;
					type[n] = s2i(steps[k]);
					k++;
					sym[n] = s2i(steps[k]);
					k++;
					a[n] = s2d(steps[k]);
					k++;
					nave[n] = s2d(steps[k]);
					k++;
					amp[n] = s2d(steps[k]);
					k++;
					namp[n] = s2d(steps[k]);
					k++;
					tx1[n] = s2d(steps[k]);
					k++;
					ty1[n] = s2d(steps[k]);
					k++;
					theta[n] = s2d(steps[k])/180.0*Math.PI;
					k++;
					xr[n] = s2d(steps[k]);
					k++;
					yr[n] = s2d(steps[k]);
					k++;
					tx2[n] = s2d(steps[k]);
					k++;
					ty2[n] = s2d(steps[k]);
					k++;
				}
				double u, v, sine, cosine, rx, ry;
				for(int j = 0; j < Ny; j++) {
					v = j*dy;
					for(int i = 0; i < Nx; i++) {
						u = i*dx;
						int n = closest(u, v, x, y, W, H);
						if(type[n] == 2) continue;
						if(type[n] < 0 || type[n] > 2) {
							System.out.println("Warning: Invalid type used for a Voronoi grain. Defaulting to overwrite.");
							type[n] = 0;
						}
						u = i*dx - tx2[n] - xr[n];
						v = j*dy - ty2[n] - yr[n];
						sine = Math.sin(-theta[n]);
						cosine = Math.cos(-theta[n]);
						rx = u*cosine - v*sine + xr[n] - tx1[n];
						ry = u*sine + v*cosine + yr[n] - ty1[n];
						arr[j][i] = type[n]*arr[j][i] + nave[n] + amp[n]*OMA(rx, ry, sym[n], a[n]) + 2.0*namp[n]*(rng.nextDouble() - 0.5);
					}
				}
			}

			// circular grain
			else if(steps[k].equals("C")) {
				k++;
				double xc = s2d(steps[k]);
				k++;
				double yc = s2d(steps[k]);
				k++;
				double R = s2d(steps[k]);
				k++;
				int type = s2i(steps[k]);
				k++;
				if(type == 2) {
					System.out.println("Warning: Phantom type used for a circular grain.");
					continue;
				}
				if(type < 0 || type > 2) {
					System.out.println("Warning: Invalid type used for a circular grain. Defaulting to overwrite.");
					type = 0;
				}
				int sym = s2i(steps[k]);
				k++;
				double a = s2d(steps[k]);
				k++;
				double nave = s2d(steps[k]);
				k++;
				double amp = s2d(steps[k]);
				k++;
				double namp = s2d(steps[k]);
				k++;
				double tx1 = s2d(steps[k]);
				k++;
				double ty1 = s2d(steps[k]);
				k++;
				double theta = s2d(steps[k])/180.0*Math.PI;
				k++;
				double xr = s2d(steps[k]);
				k++;
				double yr = s2d(steps[k]);
				k++;
				double tx2 = s2d(steps[k]);
				k++;
				double ty2 = s2d(steps[k]);
				k++;

				double u, v, x_, y_, x, y;
				double R2 = R*R;
				double sine = Math.sin(-theta);
				double cosine = Math.cos(-theta);
				int wi, wj;
				int i1 = (int)((xc - R)/dx);
				int i2 = (int)((xc + R)/dx);
				int j1 = (int)((yc - R)/dy);
				int j2 = (int)((yc + R)/dy);
				for(int j = j1; j <= j2; j++) {
					v = j*dy - yc;
					y_ = j*dy - ty2 - yr;
					wj = w(j, Ny);
					for(int i = i1; i <= i2; i++) {
						u = i*dx - xc;
						wi = w(i, Nx);
						if(u*u + v*v < R2) {
							x_ = i*dx - tx2 - xr;
							x = x_*cosine - y_*sine + xr - tx1;
							y = x_*sine + y_*cosine + yr - ty1;
							arr[wj][wi] = type*arr[wj][wi] + nave + amp*OMA(x, y, sym, a) + 2.0*namp*(rng.nextDouble() - 0.5);
						}
					}
				}
			}

			// polygon grain (only convex!)
			else if(steps[k].equals("P")) {
				k++;
				int N = s2i(steps[k]);
				k++;
				double[] xc = new double[N];
				double[] yc = new double[N];
				for(int n = 0; n < N; n++) {
					xc[n] = s2d(steps[k]);
					k++;
					yc[n] = s2d(steps[k]);
					k++;
				}
				double ax, ay, bx, by, z1, z2;
				double xmin = xc[0];
				double xmax = xc[0];
				double ymin = yc[0];
				double ymax = yc[0];
				ax = xc[0] - xc[N - 1];
				ay = yc[0] - yc[N - 1];
				bx = xc[1] - xc[0];
				by = yc[1] - yc[0];
				z1 = ax*by - ay*bx;
				for(int n = 1; n < N; n++) {
					if(xc[n] < xmin) xmin = xc[n];
					if(xc[n] > xmax) xmax = xc[n];
					if(yc[n] < ymin) ymin = yc[n];
					if(yc[n] > ymax) ymax = yc[n];
					ax = xc[n] - xc[(n - 1 + N)%N];
					ay = yc[n] - yc[(n - 1 + N)%N];
					bx = xc[(n + 1)%N] - xc[n];
					by = yc[(n + 1)%N] - yc[n];
					z2 = ax*by - ay*bx;
					if(z1*z2 < 0.0) {
						System.out.println("Error: Concave grain detected. Convex grain expected.");
						exit(1);
					}
					else z1 = z2;
				}
				int type = s2i(steps[k]);
				k++;
				if(type == 2) {
					System.out.println("Warning: Phantom type used for a polygon grain.");
					continue;
				}
				if(type < 0 || type > 2) {
					System.out.println("Warning: Invalid type used for a polygon grain. Defaulting to overwrite.");
					type = 0;
				}
				int sym = s2i(steps[k]);
				k++;
				double a = s2d(steps[k]);
				k++;
				double nave = s2d(steps[k]);
				k++;
				double amp = s2d(steps[k]);
				k++;
				double namp = s2d(steps[k]);
				k++;
				double tx1 = s2d(steps[k]);
				k++;
				double ty1 = s2d(steps[k]);
				k++;
				double theta = s2d(steps[k])/180.0*Math.PI;
				k++;
				double xr = s2d(steps[k]);
				k++;
				double yr = s2d(steps[k]);
				k++;
				double tx2 = s2d(steps[k]);
				k++;
				double ty2 = s2d(steps[k]);
				k++;

				double x, y, x_, y_, u, v;
				double sine = Math.sin(-theta);
				double cosine = Math.cos(-theta);
				int wi, wj;
				int i1 = (int)(xmin/dx);
				int i2 = (int)(xmax/dx);
				int j1 = (int)(ymin/dy);
				int j2 = (int)(ymax/dy);
				for(int j = j1; j <= j2; j++) {
					wj = w(j, Ny);
					y = j*dy;
					y_ = y - ty2 - yr;
					outer:
					for(int i = i1; i <= i2; i++) {
						wi = w(i, Nx);
						x = i*dx;
						x_ = x - tx2 - xr;
						u = x_*cosine - y_*sine + xr - tx1;
						v = x_*sine + y_*cosine + yr - ty1;

						ax = xc[0] - xc[N - 1];
						ay = yc[0] - yc[N - 1];
						bx = x - xc[N - 1];
						by = y - yc[N - 1];
						z1 = ax*by - ay*bx;
						for(int n = 1; n < N; n++) {
							ax = xc[n] - xc[n - 1];
							ay = yc[n] - yc[n - 1];
							bx = x - xc[n];
							by = y - yc[n];
							z2 = ax*by - ay*bx;
							if(z1*z2 < 0.0) continue outer;
							else z1 = z2;
						}
						arr[wj][wi] = type*arr[wj][wi] + nave + amp*OMA(u, v, sym, a) + 2.0*namp*(rng.nextDouble() - 0.5);
					}
				}

			}

			// transform
			else if(steps[k].equals("T")) {
				k++;
				double nave = s2d(steps[k]);
				k++;
				double famp = s2d(steps[k]);
				k++;
				double nave0 = 0.0;
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						nave0 += arr[j][i];
					}
				}
				nave0 *= 1.0/Nx/Ny;
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						arr[j][i] = (arr[j][i] - nave0)*famp + nave;
					}
				}
			}

			// output
			else if(steps[k].equals("O")) {
				k++;
				name = steps[k];
				k++;
				int prec = s2i(steps[k]);
				k++;
				String format = "%." + prec + "e";
				BufferedWriter writer = new BufferedWriter(new FileWriter(name));
				writer.write(Nx + " " + Ny + " ");
				writer.write(String.format(format + " " + format + "\n", dx, dy));
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						writer.write(String.format(format + "\n", arr[j][i]));
					}
				}
				writer.close();
			}

			// garbage
			else {
				System.out.println("Error: Invalid input file. Label expected instead of \"" + steps[k] + "\".");
				exit(1);
			}
		}

		reader.close();

		exit(0);

	}
}
