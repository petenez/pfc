// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package atomistic;

import java.awt.*;
import java.io.*;

import static java.lang.System.exit;

public class Colorizer {

	static double W, H;	// system dimensions

	// returns the magnitude of vector u
	public static double abs(double[] u) {
		return Math.sqrt(u[0]*u[0] + u[1]*u[1]);
	}

	// returns the argument angle of vector u
	public static double arg(double[] u) {
		return Math.atan2(u[1], u[0]);
	}

	// sums vectors u and v and stores the result in w
	public static void sum(double[] u, double[] v, double[] w) {
		w[0] = u[0] + v[0];
		w[1] = u[1] + v[1];
	}

	// subtracts vector v from u and stores the result in w
	public static void sub(double[] u, double[] v, double[] w) {
		w[0] = u[0] - v[0];
		w[1] = u[1] - v[1];
	}

	// multiplies vector u by scalar a and stores the result in w
	public static void mul(double[] u, double a, double[] w) {
		w[0] = a*u[0];
		w[1] = a*u[1];
	}

	// creates a vector of given magnitude and argument angle and stores the result in w
	public static void polar(double abs, double arg, double[] w) {
		w[0] = abs*Math.cos(arg);
		w[1] = abs*Math.sin(arg);
	}

	// wraps vector u to take periodic boundary conditions into account
	public static void wrap(double[] u) {
		if(u[0] >= 0.5*W) u[0] -= W;
		if(u[0] < -0.5*W) u[0] += W;
		if(u[1] >= 0.5*H) u[1] -= H;
		if(u[1] < -0.5*H) u[1] += H;
	}

	// saturates 0 <= d <= 1
	public static double s(double d) {
		if(d < 0.0) return 0.0;
		if(d > 1.0) return 1.0;
		return d;
	}

	// tries to parse a string to an integer
	public static int s2i(String s) {
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
	public static double s2d(String s) {
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

		// print instructions
		if(args.length == 0) {
			System.out.println("This is a tool for creating a style file where the colors of atoms are based on their bond orientations. The program expects as its arguments the input and output filenames, and the atom and bond radii. A lattice with six-fold rotational symmetry and atoms with three neighbors are expected.");
			System.out.println();
			System.out.println("An example of valid syntax:");
			System.out.println();
			System.out.println("java -jar colorizer.jar foo.xyz bar.sty 0.5 0.15");
			exit(0);
		}

		if(args.length != 4) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

		String input = args[0];		// input filename
		String output = args[1];	// output filename
		double R = s2d(args[2]);	// atom radius
		double r = s2d(args[3]);	// bond radius

		if(!new File(input).exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}

		int N;	// # of atoms
		BufferedReader reader = new BufferedReader(new FileReader(input));
		String line = reader.readLine();
		String[] words = line.trim().split("\\s+");
		W = s2d(words[0]);
		H = s2d(words[1]);
		for(N = 0; reader.readLine() != null; N++);	// count # of atoms
		reader.close();

		double[][] xyz = new double[N][3];	// atom coordinates
		int[][] nl = new int[N][3];			// neighbor lists

		reader = new BufferedReader(new FileReader(input));
		reader.readLine();
		for(int n = 0; n < N; n++) {
			words = reader.readLine().trim().split("\\s+");
			xyz[n][0] = s2d(words[1]);
			xyz[n][1] = s2d(words[2]);
			xyz[n][2] = s2d(words[3]);
			nl[n][0] = s2i(words[4]);
			nl[n][1] = s2i(words[5]);
			nl[n][2] = s2i(words[6]);
		}
		reader.close();

		double _3 = 1.0/3.0;
		double pi_3 = Math.PI/3.0;
		double pi_6 = Math.PI/6.0;
		double _pi_3 = 3.0/Math.PI;
		double _255 = 1.0/255.0;
		double arg, ARG, abs, red, green, blue;
		double[] d = new double[2];		// difference vector from nth atom to its mth neighbor
		double[] p = new double[2];		// polar coordinates vector
		double[] D;						// sum of modified d vectors
		Color c;						// color object for HSB-to-RGB color transformation
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		for(int n = 0; n < N; n++) {
			D = new double[2];								// sum zeroed
			for(int m = 0; m < 3; m++) {
				sub(xyz[nl[n][m]], xyz[n], d);				// difference vector
				wrap(d);									// periodic boundary conditions
				arg = arg(d);								// argument angle
				if(abs(D) == 0.0) {							// first d to be added to the sum?
					while(arg >= pi_3) arg -= pi_3;			// rotate 0 <= arg <= pi/3 (six-fold symmetry assumed)
					while(arg < 0.0) arg += pi_3;
				}
				else {
					ARG = arg(D);
					while(arg - ARG >= pi_6) arg -= pi_3;	// rotate so that |arg - ARG| <= pi/6
					while(arg - ARG < -pi_6) arg += pi_3;
				}
				polar(1.0, arg, p);							// unit vector with new argument angle
				sum(D, p, D);								// add to sum
			}
			mul(D, _3, D);									// normalize
			abs = s(abs(D));								// magnitude 0 <= abs <= 1 (in case of floating-point errors)
			arg = arg(D)*_pi_3;								// normalize
			if(arg >= 1.0) arg = arg%1.0;					// 0 <= arg <= 1
			if(arg < 0.0) arg += 1.0;
			arg = (arg + 0.5)%1.0;							// shift to match colors with those of orientation fields
			c = new Color(Color.HSBtoRGB((float)arg, (float)abs, 1));	// new color by hue, saturation and brightness
			red = s(c.getRed()*_255);						// 0 <= red < 256 --> 0 <= red <= 1
			green = s(c.getGreen()*_255);					// ...
			blue = s(c.getBlue()*_255);
			writer.write(R + " " + r + " " + red + " " + green + " " + blue + " 0.0 \n");
		}
		writer.close();

	}
}
