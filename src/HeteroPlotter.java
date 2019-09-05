// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package visualization;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import static java.lang.System.exit;

public class HeteroPlotter {

	// tries to parse a string to an integer
	static int s2i(String s) {
		int i = 0;
		try {
			i = Integer.parseInt(s);
		} catch (NumberFormatException e) {
			System.out.println("Error: Invalid String. Integer expected instead of \"" + s + "\".");
			exit(1);
		}
		return i;
	}

	// tries to parse a string to a double
	static double s2d(String s) {
		double d = 0.0;
		try {
			d = Double.parseDouble(s);
		} catch (NumberFormatException e) {
			System.out.println("Error: Invalid String. Decimal number expected instead of \"" + s + "\".");
			exit(1);
		}
		return d;
	}

	public static void main(String[] args) throws IOException {

		int A = args.length;
		if(A == 0) {
			System.out.println("This is a tool for visualizing PFC heterostructures (see arXiv:1908.05564). As its arguments, the program expects equal numbers of density field data files and corresponding smoothed density field data files plus an output data file. Note that at least two density fields and two corresponding smoothed density fields are expected.");
			System.out.println();
			System.out.println("An example of valid syntax:");
			System.out.println();
			System.out.println("java -jar heteroplotter.jar foo.n0 foo.n1 bar.eta0 bar.eta1 result.png");
			exit(0);
		}
		else if(A%2 != 1 || A < 5) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

		int N = A/2;
		for(int n = 0; n < N; n++) {
			if(!new File(args[n]).exists() || !new File(args[N + n]).exists()) {
				System.out.println("Error: Input data file not found.");
				exit(1);
			}
		}

		int Nx, Ny;
		BufferedReader[] nreaders = new BufferedReader[N];
		BufferedReader[] etareaders = new BufferedReader[N];
		nreaders[0] = new BufferedReader(new FileReader(args[0]));
		String[] words = nreaders[0].readLine().trim().split("\\s+");
		nreaders[0].close();
		Nx = s2i(words[0]);		// dimensions of first data file
		Ny = s2i(words[1]);
		// compare with the dimensions of the rest
		for(int n = 0; n < N; n++) {
			nreaders[n] = new BufferedReader(new FileReader(args[n]));
			words = nreaders[n].readLine().trim().split("\\s+");
			if(s2i(words[0]) != Nx || s2i(words[1]) != Ny) {
				System.out.println("Error: Mismatching data dimensions.");
				exit(1);
			}
			etareaders[n] = new BufferedReader(new FileReader(args[N + n]));
			words = etareaders[n].readLine().trim().split("\\s+");
			if(s2i(words[0]) != Nx || s2i(words[1]) != Ny) {
				System.out.println("Error: Mismatching data dimensions.");
				exit(1);
			}
		}

		double[][][] ns = new double[Ny][Nx][N];	// array for density data
		double[][][] etas = new double[Ny][Nx][N];	// array for smoothed density data
		double[] naves = new double[N];				// average densities
		double[] nmins = new double[N];				// density minima
		Arrays.fill(nmins, 1.0e300);
		double[] nmaxs = new double[N];				// density maxima
		Arrays.fill(nmaxs, -1.0e300);
		double[] etamins = new double[N];			// smoothed density minima
		Arrays.fill(etamins, 1.0e300);
		double[] etamaxs = new double[N];			// smoothed density maxima
		Arrays.fill(etamaxs, -1.0e300);
		// read data
		for(int j = 0; j < Ny; j++) {
			for(int i = 0; i < Nx; i++) {
				for(int n = 0; n < N; n++) {
					etas[j][i][n] = Double.parseDouble(etareaders[n].readLine());
					if(etas[j][i][n] < etamins[n]) etamins[n] = etas[j][i][n];
					if(etas[j][i][n] > etamaxs[n]) etamaxs[n] = etas[j][i][n];
					ns[j][i][n] = Double.parseDouble(nreaders[n].readLine());
					naves[n] += ns[j][i][n];
					ns[j][i][n] -= etas[j][i][n];
					if(ns[j][i][n] < nmins[n]) nmins[n] = ns[j][i][n];
					if(ns[j][i][n] > nmaxs[n]) nmaxs[n] = ns[j][i][n];
				}
			}
		}
		for(int n = 0; n < N; n++) naves[n] *= 1.0/Nx/Ny;
		double abs;		// magnitude
		double arg;		// phase
		double mask;	// mask
		double twopi_N = 2.0*Math.PI/N;
		double absmax = -1.0e300;
		double mmax = -1.0e300;
		double[][] m = new double[Ny][Nx];
		double[][][] v = new double[Ny][Nx][2];
		for(int j = 0; j < Ny; j++) {
			for(int i = 0; i < Nx; i++) {
				for(int n = 0; n < N; n++) {
					// mask is zero for a constant density and unity for a crystalline density
					if(naves[n] < 0.0) mask = (etas[j][i][n] - etamins[n])/(etamaxs[n] - etamins[n]);
					else mask = (etas[j][i][n] - etamaxs[n])/(etamins[n] - etamaxs[n]);
					// density is first mapped between zero and unity and then masked
					abs = mask*(ns[j][i][n] - nmins[n])/(nmaxs[n] - nmins[n]);
					// each density field is given a phase (corresponding to a hue)
					arg = twopi_N*n;
					// combined hue-saturation vector
					v[j][i][0] += abs*Math.cos(arg);
					v[j][i][1] += abs*Math.sin(arg);
					// combined density magnitude
					m[j][i] += abs*abs;
				}
				// magnitude of vector
				abs = Math.sqrt(v[j][i][0]*v[j][i][0] + v[j][i][1]*v[j][i][1]);
				// magnitude (Pythagorean)
				m[j][i] = Math.sqrt(m[j][i]);
				if(abs > absmax) absmax = abs;
				if(m[j][i] > mmax) mmax = m[j][i];
			}
		}
		double _absmax = 1.0/absmax;
		double _mmax = 1.0/mmax;
		double _2pi = 0.5/Math.PI;
		float h, s, b;	// hue, saturation, brightness
		BufferedImage image = new BufferedImage(Nx, Ny, BufferedImage.TYPE_INT_RGB);
		for(int j = 0; j < Ny; j++) {
			for(int i = 0; i < Nx; i++) {
				// phase of hue-saturation vector --> hue
				h = (float)(Math.atan2(v[j][i][1], v[j][i][0])*_2pi + 0.5);
				// magnitude of hue-saturation vector --> saturation
				s = (float)(Math.sqrt(v[j][i][0]*v[j][i][0] + v[j][i][1]*v[j][i][1])*_absmax);
				// density magnitude --> brightness
				b = (float)(m[j][i]*_mmax);
				image.setRGB(i, Ny - j - 1, Color.HSBtoRGB(h, s, b));
			}
		}
		ImageIO.write(image, "png", new File(args[A - 1]));
	}
}
