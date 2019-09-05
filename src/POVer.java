// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package atomistic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import static java.lang.System.exit;

// a simple container for a light source
class Light {
	double[] r = new double[3];	// location vector
	double[] c = new double[3];	// RGB color vector

	Light(double x, double y, double z, double red, double green, double blue) {
		r[0] = x;
		r[1] = y;
		r[2] = z;
		c[0] = red;
		c[1] = green;
		c[2] = blue;
	}
}

public class POVer {

	// saturates 0 <= d <= 1
	private static double s(double d) {
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
			System.out.println("Error: Invalid parameter file. Integer expected instead of \"" + s + "\".");
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
			System.out.println("Error: Invalid parameter file. Decimal number expected instead of \"" + s + "\".");
			exit(1);
		}
		return d;
	}

	// returns the magnitude of a vector
	public static double abs(double[] a) {
		return Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	}

	// returns the division of a vector
	public static double[] div(double[] a, double b) {
		return new double[] {a[0]/b, a[1]/b, a[2]/b};
	}

	// returns the cross product of two vectors
	public static double[] x(double[] a, double[] b) {
		return new double[] {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
	}

	public static void main(String[] args) throws IOException {

		if(args.length == 0) {
			System.out.println("This is a tool for generating POV-Ray scripts for rendering atomic configurations from an xyz-file. The program expects one argument specifying a parameter file or, alternatively, three arguments: input and output filenames and an atom radius. The latter results in all atoms being black and having the same radius. The view angle will be zero, corresponding to an orthographic view.\n" +
									   "\n" +
									   "An example of a simple parameter file:\n" +
									   "\n" +
									   "# input file\n" +
									   "I foo.xyz\n" +
									   "\n" +
									   "# output file\n" +
									   "O bar.pov\n" +
									   "\n" +
									   "# camera position (x y z)\n" +
									   "C 0.0 -0.001 1000.0\n" +
									   "\n" +
									   "# \"look at\" point (x y z)\n" +
									   "P 0.0 0.0 0.0\n" +
									   "\n" +
									   "# light sources (x y z red green blue)\n" +
									   "L 0.0 0.0 10000.0 1.0 1.0 1.0\n" +
									   "L 0.0 10000.0 0.0 1.0 0.0 0.0\n" +
									   "\n" +
									   "# background color (red green blue transparency)\n" +
									   "B 1.0 1.0 1.0 0.0\n" +
									   "\n" +
									   "# view angle (degrees)\n" +
									   "A 30.0\n" +
									   "\n" +
									   "# styles for atom types (type atom_radius bond_radius red green blue transparency)\n" +
									   "S C 0.5 0.15 0.0 1.0 1.0 0.0\n" +
									   "S B 0.5 0.15 1.0 0.0 1.0 0.0\n" +
									   "S N 0.5 0.15 1.0 1.0 0.0 0.0\n" +
									   "\n" +
									   "# note that \"*\" can be used to specify the style for all unspecified atom types\n" +
									   "S * 0.5 0.15 0.0 0.0 0.0 0.5\n" +
									   "\n" +
									   "# style file (style for each atom separately; if uncommented, this will override the specifications above)\n" +
									   "#F test.sty\n" +
									   "\n" +
									   "This example demonstrates two light sources - one white, one red - and a white, non-transparent background. The \"C\", \"B\" and \"N\" atom types are specified as cyan, magenta and yellow, and all others as semi-transparent black.\n" +
									   "\n" +
									   "Warning: I never fully understood the logic behind POV-Ray's coordinate system and camera vectors. Try to avoid having the camera directly above the \"look at\" point, for example.");
			exit(0);
		}

		String param = null;	// file name for parameter file
		String input = null;	// file name for input file
		String style = null;	// file name for style file
		String output = null;	// file name for output file
		BufferedReader reader;

		double W;	// system width
		double H;	// system height
		int N;		// # of atoms

		String[] types;		// atom types
		double[][] xyz;		// atom coordinates
		int[][] nl;			// neighbor lists
		double[] R;			// atom radii
		double[] r;			// bond radii
		double[][] rgbt;	// atom colors (red, green, blue, transparency)
		// a hashmap connecting an atom's type to corresponding parameters (radii and color)
		HashMap<String, double[]> styles = new HashMap<String, double[]>();

		double alpha = 0.0;	// view angle
		double[] c = null;	// camera location vector
		double[] p = null;	// "look at" point
		double[] b = null;	// background color
		ArrayList<Light> lights = new ArrayList<Light>(1);	// list of light sources

		boolean b_quick = false;	// quick run?

		// load parameters from a file
		if(args.length == 1) {
			param = args[0];
			File file = new File(param);
			if(!file.exists()) {
				System.out.println("Error: Parameter file not found.");
				exit(1);
			}
			String line;	// line of input
			String[] words;	// line split into "words" delimited by any whitespace characters
			ArrayList<String> steps = new ArrayList<String>();
			reader = new BufferedReader(new FileReader(param));

			// clean input
			outer:
			while((line = reader.readLine()) != null) {
				words = line.trim().split("\\s+");
				for(int w = 0; w < words.length; w++) {
					if(words[w].isEmpty()) continue;				// empty
					if(words[w].charAt(0) == '#') continue outer;	// comment
					steps.add(words[w]);
				}
			}

			int k = 0;
			while(k < steps.size()) {	// actually load parameters

				// set camera position
				if(steps.get(k).equals("C")) {
					k++;						// move to next "word" of parameter file
					c = new double[3];
					c[0] = s2d(steps.get(k));
					k++;
					c[1] = s2d(steps.get(k));
					k++;
					c[2] = s2d(steps.get(k));
					k++;
				}

				// set "look at" point
				if(steps.get(k).equals("P")) {
					k++;
					p = new double[3];
					p[0] = s2d(steps.get(k));
					k++;
					p[1] = s2d(steps.get(k));
					k++;
					p[2] = s2d(steps.get(k));
					k++;
				}

				// add light source
				if(steps.get(k).equals("L")) {
					k++;
					double x = s2d(steps.get(k));
					k++;
					double y = s2d(steps.get(k));
					k++;
					double z = s2d(steps.get(k));
					k++;
					double red = s2d(steps.get(k));
					k++;
					double green = s2d(steps.get(k));
					k++;
					double blue = s2d(steps.get(k));
					k++;
					lights.add(new Light(x, y, z, red, green, blue));
				}

				// set background color
				else if(steps.get(k).equals("B")) {
					k++;
					b = new double[4];
					b[0] = s2d(steps.get(k));
					k++;
					b[1] = s2d(steps.get(k));
					k++;
					b[2] = s2d(steps.get(k));
					k++;
					b[3] = s2d(steps.get(k));	// can be transparent
					k++;
				}

				// set view angle
				else if(steps.get(k).equals("A")) {
					k++;
					alpha = s2d(steps.get(k));
					k++;
					if(alpha < 0.0) {
						System.out.println("Warning: Negative view angle specified. A value from the range [0, 180) degrees is expected. Setting the view angle to zero.");
						alpha = 0.0;
					}
					else if(alpha >= 180.0) {
						System.out.println("Warning: View angle >= 180 degrees specified. A value from the range [0, 180) degrees is expected. Setting the view angle to zero.");
						alpha = 0.0;
					}
				}

				// specify input file
				else if(steps.get(k).equals("I")) {
					k++;
					input = steps.get(k);
					k++;
				}

				// specify output file
				else if(steps.get(k).equals("O")) {
					k++;
					output = steps.get(k);
					k++;
				}

				// specify a style file specifying the style of each atom separately
				else if(steps.get(k).equals("F")) {
					k++;
					style = steps.get(k);
					k++;
				}

				// specify styles for atom types
				// note that a style can be applied to all unspecified atom types by using "*":
				// S * 0.5 0.15 0.0 0.0 0.0 0.0
				else if(steps.get(k).equals("S")) {
					k++;
					String type = steps.get(k);
					k++;
					double[] arr = new double[6];
					arr[0] = s2d(steps.get(k));		// atom radius
					k++;
					if(arr[0] < 0.0) arr[0] = 0.0;
					arr[1] = s2d(steps.get(k));		// bond radius
					k++;
					if(arr[1] < 0.0) arr[1] = 0.0;
					arr[2] = s(s2d(steps.get(k)));	// red
					k++;
					arr[3] = s(s2d(steps.get(k)));	// green
					k++;
					arr[4] = s(s2d(steps.get(k)));	// blue
					k++;
					arr[5] = s(s2d(steps.get(k)));	// transparency
					k++;
					styles.put(type, arr);
				}

				// garbage
				else {
					System.out.println("Error: Invalid parameter file. Label expected instead of \"" + steps.get(k) + "\".");
					exit(1);
				}
			}
		}
		// quick run
		else if(args.length == 3) {
			// specify input and output files and atom radius
			input = args[0];
			output = args[1];
			double R_ = s2d(args[2]);
			if(R_ < 0.0) {
				System.out.println("Warning: Negative atom radius specified. Positive value expected. Radius set to zero.");
				R_ = 0.0;
			}
			double[] arr = new double[] {R_, 0.0, 0.0, 0.0, 0.0, 0.0};	// no bonds, black atoms
			styles.put("*", arr);	// style for all atoms
			b_quick = true;
		}
		// garbage
		else {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

		// input file exists?
		File file = new File(input);
		if(!file.exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}

		// extract W, H and N (# of atoms)
		reader = new BufferedReader(new FileReader(input));
		String line = reader.readLine();
		String[] words = line.trim().split("\\s+");
		W = s2d(words[0]);
		H = s2d(words[1]);
		for(N = 0; reader.readLine() != null; N++);
		reader.close();

		// initialize arrays
		types = new String[N];
		xyz = new double[N][3];
		nl = new int[N][];
		R = new double[N];
		r = new double[N];
		rgbt = new double[N][4];

		// read types, coordinates and neighbors
		reader = new BufferedReader(new FileReader(input));
		reader.readLine();
		int M;
		for(int n = 0; n < N; n++) {
			words = reader.readLine().trim().split("\\s+");
			types[n] = words[0];
			xyz[n][0] = s2d(words[1]);
			xyz[n][1] = s2d(words[2]);
			xyz[n][2] = s2d(words[3]);
			M = words.length - 4;	// # of neighbors
			nl[n] = new int[M];		// initialize nth neighbor list
			for(int m = 0; m < M; m++) {
				nl[n][m] = s2i(words[m + 4]);
			}
		}
		reader.close();

		// read style file
		// file name specified and corresponding style file exists?
		if(style != null && new File(style).exists()) {
			reader = new BufferedReader(new FileReader(style));
			for(int n = 0; n < N; n++) {
				words = reader.readLine().trim().split("\\s+");
				R[n] = s2d(words[0]);
				if(R[n] < 0.0) {
					System.out.println("Warning: Negative atom radius specified. Positive value expected. Radius set to zero.");
					R[n] = 0.0;
				}
				r[n] = s2d(words[1]);
				if(r[n] < 0.0) {
					System.out.println("Warning: Negative bond radius specified. Positive value expected. Radius set to zero.");
					r[n] = 0.0;
				}
				rgbt[n][0] = s(s2d(words[2]));
				rgbt[n][1] = s(s2d(words[3]));
				rgbt[n][2] = s(s2d(words[4]));
				rgbt[n][3] = s(s2d(words[5]));
			}
			reader.close();
		}
		// style file specified but not found
		else if(style != null) {
			System.out.println("Error: Style file not found.");
			exit(1);
		}
		// quick run or atom styles specified in parameter file
		else if(!styles.isEmpty()) {
			double[] arr;
			for(int n = 0; n < N; n++) {
				arr = styles.get(types[n]);
				// style for given type not specified? --> use type "*"
				if(arr == null) arr = styles.get("*");
				// type "*" not specified?
				if(arr == null) {
					System.out.println("Warning: Unknown atom type \"" + types[n] + "\" detected. Set invisible.");
					arr = new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
				}
				R[n] = arr[0];
				r[n] = arr[1];
				rgbt[n][0] = arr[2];
				rgbt[n][1] = arr[3];
				rgbt[n][2] = arr[4];
				rgbt[n][3] = arr[5];
			}
		}
		else {
			System.out.println("Error: No styles specified.");
			exit(1);
		}

		// camera/"look at"/background/light sources not specified
		if(c == null) {
			if(!b_quick) {
				System.out.println("Warning: Camera location not specified. Using default location and setting view angle to zero.");
				alpha = 0.0;	// --> orthographic view
			}
			c = new double[] {0.5*W, 0.5*H, Math.max(W, H)};
		}
		if(p == null) {
			if(!b_quick) System.out.println("Warning: Look at point not specified. Using default location.");
			p = new double[] {0.5*W, 0.5*H, 0.0};
		}
		if(b == null) {
			if(!b_quick) System.out.println("Warning: Background color not specified. Using default color.");
			b = new double[] {1.0, 1.0, 1.0, 0.0};
		}
		if(lights.size() == 0) {
			if(!b_quick) System.out.println("Warning: No light sources specified. Adding default light source.");
			double x = p[0] + 1.0e6*(c[0] - p[0]);
			double y = p[1] + 1.0e6*(c[1] - p[1]);
			double z = p[2] + 1.0e6*(c[2] - p[2]);
			lights.add(new Light(x, y, z, 1.0, 1.0, 1.0));
		}

		// write pov file
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		writer.write(
				"#version 3.7;\n" +
						"global_settings{assumed_gamma 1.0}\n" +
						//"#default{finish{ambient 0.1 diffuse 0.9}}\n" +
						"\n" +
						"#include \"colors.inc\"\n\n");

		// specify camera
		double sign = 1.0;	// a hack (I don't understand the full logic behind POV-Ray's coordinate system and camera vectors, but this seems to make things work)
		// note also: left-handed to right-handed coordinates --> y and z swapped
		if(alpha == 0.0) {
			double dx = c[0] - p[0];
			double dy = c[1] - p[1];
			double dz = c[2] - p[2];
			double abs = Math.sqrt(dx*dx + dy*dy + dz*dz);
			// view angle set to fit the larger dimension to the picture
			alpha = Math.atan(0.5*Math.max(W, H)/abs)/Math.PI*360.0;
			writer.write(
					"camera{ orthographic angle " + alpha + "\n" +
							"	location < " + c[0] + ", " + c[2] + ", " + sign*c[1]  + " >\n" +
							"	look_at < " + p[0] + ", " + p[2] + ", " + sign*p[1] + " >\n" +
							"	right image_width/image_height*x\n" +
							"}\n\n"
			);
		}
		else {
			double[] k = new double[] {0.0, 0.0, 1.0};	// unit vector in z direction
			// vector from "look at" point to camera location
			double[] p_c = new double[] {c[0] - p[0], c[1] - p[1], c[2] - p[2]};
			if(p_c[1] >= 0.0) sign = -1.0;
			double[] x = x(k, p_c);	// cross product of k and p_c (horizontal direction in picture always parallel to the x-y plane)
			double absx = abs(x);
			if(absx == 0.0) x = new double[] {1.0, 0.0, 0.0};	// looking from above? --> x to the direction of the positive x axis
			else x = div(x, absx);	// make unit vector
			writer.write(
					"camera{ perspective angle " + alpha + "\n" +
							"	location < " + c[0] + ", " + c[2] + ", " + sign*c[1]  + " >\n" +
							"	look_at < " + p[0] + ", " + p[2] + ", " + sign*p[1] + " >\n" +
							"	right image_width/image_height*< " + x[0] + ", " + x[2] + ", " + sign*x[1] + " >\n" +
							"}\n\n"
			);
		}

		// specify light sources
		Light l;
		for(int n = 0; n < lights.size(); n++) {
			l = lights.get(n);
			writer.write(
					"light_source{< " + l.r[0] + ", " + l.r[2] + ", " + sign*l.r[1] + " > color srgb < " + l.c[0] + " " + l.c[1] + " " + l.c[2] + " >}\n\n"
			);
		}

		// specify background color
		writer.write("background{color srgbt < " + b[0] + " " + b[1] + " " + b[2] + " " + b[3] + " >}\n\n");

		// specify atoms and bonds
		double a = 0.25;	// coefficient for ambient light (a little, not too much)
		double x1, y1, z1, x2, y2, z2, xa, ya, za, ra;
		double[] rgbt_ = new double[4];
		for(int n = 0; n < N; n++) {					// nth atom
			x1 = xyz[n][0];								// coordinates
			y1 = xyz[n][1];
			z1 = xyz[n][2];
			for(int m = 0; m < nl[n].length; m++) {		// mth neighbor
				x2 = xyz[nl[n][m]][0];					// coordinates
				y2 = xyz[nl[n][m]][1];
				z2 = xyz[nl[n][m]][2];
				if(x2 - x1 < -0.5*W) x2 += W;			// periodic boundaries (bond stubs)
				if(x2 - x1 > 0.5*W) x2 -= W;
				if(y2 - y1 < -0.5*H) y2 += H;
				if(y2 - y1 > 0.5*H) y2 -= H;
				xa = 0.5*(x1 + x2);						// average coordinates (bond midpoints)
				ya = 0.5*(y1 + y2);
				za = 0.5*(z1 + z2);
				ra = 0.5*(r[n] + r[nl[n][m]]);			// average bond thickness
				rgbt_[0] = 0.5*(rgbt[n][0] + rgbt[nl[n][m]][0]);	// average colors
				rgbt_[1] = 0.5*(rgbt[n][1] + rgbt[nl[n][m]][1]);
				rgbt_[2] = 0.5*(rgbt[n][2] + rgbt[nl[n][m]][2]);
				rgbt_[3] = 0.5*(rgbt[n][3] + rgbt[nl[n][m]][3]);
				writer.write(	// bonds
						"cone{\n" +
								"	< " + x1 + ", " + z1 + ", " + sign*y1 + " >, " + r[n] + "\n" +
								"	< " + xa + ", " + za + ", " + sign*ya + " >, " + ra + "\n" +
								"	open\n" +
								"	texture{\n" +
								"		pigment{color srgbt < " + rgbt_[0] + " " + rgbt_[1] + " " + rgbt_[2] + " " + rgbt_[3] + " >}\n" +
								"		finish{ambient srgb < " + a*b[0] + " " + a*b[1] + " " + a*b[2] + " > specular 0.1 phong 0.25}\n" +
								"	}\n" +
								"}\n\n");
			}
			writer.write(		// atoms
					"sphere{\n" +
							"	< " + x1 + ", " + z1 + " , " + sign*y1 + " >, " + R[n] + "\n" +
							"	texture{\n" +
							"		pigment{color srgbt < " + rgbt[n][0] + " " + rgbt[n][1] + " " + rgbt[n][2] + " " + rgbt[n][3] + " >}\n" +
							"		finish{ambient srgb < " + a*b[0] + " " + a*b[1] + " " + a*b[2] + " > specular 0.1 phong 0.25}\n" +
							"	}\n" +
							"}\n\n");
		}
		reader.close();

		writer.close();
	}
}