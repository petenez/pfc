// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package data;

import java.io.*;
import static java.lang.System.exit;

public class Manipulator {

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

	// wraps negative indices, checks they are within bounds
	static void indices(int[] i) {
		if(i[0] < 0) i[0] += i[2];
		if(i[0] >= i[2]) {
			System.out.println("Error: Index out of bounds.");
			exit(1);
		}
		if(i[1] < 0) i[1] += i[2];
		if(i[1] >= i[2]) {
			System.out.println("Error: Index out of bounds.");
			exit(1);
		}
		if(i[0] > i[1]) {
			System.out.println("Error: Start index > end index.");
		}
	}

	// extracts x dimension of data from file "name"
	static int Nx(String name) throws IOException {
		if(!new File(name).exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}
		BufferedReader reader = new BufferedReader(new FileReader(name));
		int Nx = s2i(reader.readLine().trim().split("\\s+")[0]);
		reader.close();
		return Nx;
	}

	// extracts y dimension of data from file "name"
	static int Ny(String name) throws IOException {
		if(!new File(name).exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}
		BufferedReader reader = new BufferedReader(new FileReader(name));
		int Ny = s2i(reader.readLine().trim().split("\\s+")[1]);
		reader.close();
		return Ny;
	}

	// extracts x discretization of data from file "name"
	static double dx(String name) throws IOException {
		if(!new File(name).exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}
		BufferedReader reader = new BufferedReader(new FileReader(name));
		double dx = s2d(reader.readLine().trim().split("\\s+")[2]);
		reader.close();
		return dx;
	}

	// extracts y discretization of data from file "name"
	static double dy(String name) throws IOException {
		if(!new File(name).exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}
		BufferedReader reader = new BufferedReader(new FileReader(name));
		double dy = s2d(reader.readLine().trim().split("\\s+")[3]);
		reader.close();
		return dy;
	}

	// extracts # of data columns from file "name"
	static int cols(String name) throws IOException {
		if(!new File(name).exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}
		BufferedReader reader = new BufferedReader(new FileReader(name));
		reader.readLine();
		int cols = reader.readLine().trim().split("\\s+").length;
		reader.close();
		return cols;
	}

	// shorthand for initializing a BufferedReader
	static BufferedReader read(String name) throws IOException {
		return new BufferedReader(new FileReader(name));
	}

	// shorthand for initializing a BufferedWriter
	static BufferedWriter write(String name) throws IOException {
		return new BufferedWriter(new FileWriter(name));
	}

	public static void main(String[] args) throws IOException {

		if(args.length == 0) {
			System.out.println("This is a tool for manipulating PFC and other similar data files. It allows creating constant-value data, extracting and inserting data (overwriting or appending or both), rotating and flipping data, and applying various mathematical operations to data.");
			System.out.println();
			System.out.println("Examples of valid syntax:");
			System.out.println();
			System.out.println("# constant");
			System.out.println("# optional arguments are in brackets");
			System.out.println("java -jar manipulator.jar c output.n 100 100 1 [1.0]");
			System.out.println();
			System.out.println("# extraction");
			System.out.println("# default to extract everything; inclusive end indices; negative values wrapped");
			System.out.println("java -jar manipulator.jar e input.n output.n [i 0 10] [j 10 20] [k 0 -1]");
			System.out.println();
			System.out.println("# insertion");
			System.out.println("# default start indices 0 0 0");
			System.out.println("java -jar manipulator.jar i from.n to.n output.n [1 2 3]");
			System.out.println();
			System.out.println("# rotation");
			System.out.println("# valid rotations: +-90, +-180, +-270 degrees; default 90 degrees");
			System.out.println("java -jar manipulator.jar r input.n output.n [90]");
			System.out.println();
			System.out.println("# flip");
			System.out.println("# valid flip directions: x and h, and y and v");
			System.out.println("java -jar manipulator.jar f input.n output.n [x]");
			System.out.println();
			System.out.println("# addition");
			System.out.println("# alternatives indicated with slash");
			System.out.println("java -jar manipulator.jar add input.n 1.0/another.n output.n");
			System.out.println();
			System.out.println("# subtraction");
			System.out.println("java -jar manipulator.jar sub input.n 1.0/another.n output.n");
			System.out.println();
			System.out.println("# multiplication");
			System.out.println("java -jar manipulator.jar mul input.n 1.0/another.n output.n");
			System.out.println();
			System.out.println("# division");
			System.out.println("java -jar manipulator.jar div input.n 1.0/another.n output.n");
			System.out.println();
			System.out.println("# minimum");
			System.out.println("java -jar manipulator.jar min input.n 1.0/another.n output.n");
			System.out.println();
			System.out.println("# maximum");
			System.out.println("java -jar manipulator.jar max input.n 1.0/another.n output.n");
			System.out.println();
			System.out.println("# modulo");
			System.out.println("java -jar manipulator.jar mod input.n 3.14 output.n");
			System.out.println();
			System.out.println("# power");
			System.out.println("java -jar manipulator.jar pow input.n 2.0 output.n");
			System.out.println();
			System.out.println("# exponential");
			System.out.println("# base can be specified; default e");
			System.out.println("java -jar manipulator.jar exp input.n [10.0] output.n");
			System.out.println();
			System.out.println("# logarithm");
			System.out.println("# base can be specified; default e");
			System.out.println("java -jar manipulator.jar log input.n [2.0] output.n");
			System.out.println();
			System.out.println("# sine");
			System.out.println("java -jar manipulator.jar sin input.n output.n");
			System.out.println();
			System.out.println("# cosine");
			System.out.println("java -jar manipulator.jar cos input.n output.n");
			System.out.println();
			System.out.println("# tangent");
			System.out.println("java -jar manipulator.jar tan input.n output.n");
			System.out.println();
			System.out.println("# arcsine");
			System.out.println("java -jar manipulator.jar asin input.n output.n");
			System.out.println();
			System.out.println("# arccosine");
			System.out.println("java -jar manipulator.jar acos input.n output.n");
			System.out.println();
			System.out.println("# arctangent");
			System.out.println("java -jar manipulator.jar atan input.n output.n");
			System.out.println();
			System.out.println("# hyperbolic sine");
			System.out.println("java -jar manipulator.jar sinh input.n output.n");
			System.out.println();
			System.out.println("# hyperbolic cosine");
			System.out.println("java -jar manipulator.jar cosh input.n output.n");
			System.out.println();
			System.out.println("# hyperbolic tangent");
			System.out.println("java -jar manipulator.jar tanh input.n output.n");
			System.out.println();
			System.out.println("# complex conjugate");
			System.out.println("java -jar manipulator.jar conj input.n output.n");
			System.out.println();
			System.out.println("# argument");
			System.out.println("java -jar manipulator.jar arg input.n output.n");
			System.out.println();
			System.out.println("# magnitude");
			System.out.println("java -jar manipulator.jar abs input.n output.n");

			exit(0);
		}

		int s = 0;

		// constant
		// creates data with a constant value and given dimensions
		if(args[0].toLowerCase().equals("c")) {
			s++;

			String output = args[s++];	// name for output data file
			int Nx = s2i(args[s++]);	// data dimensions
			int Ny = s2i(args[s++]);
			int cols = s2i(args[s++]);	// data columns
			double d = 0.0;
			if(s < args.length) d = s2d(args[s++]);	// zero unless specified otherwise

			BufferedWriter writer = write(output);
			writer.write(Nx + " " + Ny + " 1.0 1.0\n");	// header
			for(int j = 0; j < Ny; j++) {
				for(int i = 0; i < Nx; i++) {
					for(int k = 0; k < cols; k++) {
						writer.write(d + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// extract
		// extracts data from a data file
		else if(args[0].toLowerCase().equals("e")) {
			s++;
			String input = args[s++];	// name for input data file
			String output = args[s++];

			// extraction bounds (inclusive end indices)
			// default is to extract everything
			int i0 = 0;					// start x coordinate
			int i1 = -1;				// end x coordinate (negative ones wrapped)
			int j0 = 0;					// y
			int j1 = -1;
			int k0 = 0;					// column
			int k1 = -1;

			while(s < args.length) {	// bounds specified?
				switch(args[s++]) {
					case "i":
						i0 = s2i(args[s++]);
						i1 = s2i(args[s++]);
						break;
					case "j":
						j0 = s2i(args[s++]);
						j1 = s2i(args[s++]);
						break;
					case "k":
						k0 = s2i(args[s++]);
						k1 = s2i(args[s++]);
						break;
				}
			}

			// dimensions, discretizations and # of columns extracted
			int Nx = Nx(input);
			int Ny = Ny(input);
			double dx = dx(input);
			double dy = dy(input);
			int cols = cols(input);

			// bounds processed
			int[] temp;
			indices(temp = new int[] {i0, i1, Nx});
			i0 = temp[0];
			i1 = temp[1];
			indices(temp = new int[] {j0, j1, Ny});
			j0 = temp[0];
			j1 = temp[1];
			indices(temp = new int[] {k0, k1, cols});
			k0 = temp[0];
			k1 = temp[1];

			String line;
			String[] words;
			// array for extracted data
			double[][][] data = new double[j1 - j0 + 1][i1 - i0 + 1][k1 - k0 + 1];

			BufferedReader reader = read(input);
			reader.readLine();	// skip header
			for(int j = 0; j < Ny; j++) {
				for(int i = 0; i < Nx; i++) {
					line = reader.readLine();
					if(i < i0 || i > i1 || j < j0 || j > j1) continue;	// out of bounds?
					words = line.trim().split("\\s+");
					for(int k = k0; k <= k1; k++) {
						data[j - j0][i - i0][k - k0] = s2d(words[k]);
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write((i1 - i0 + 1) + " " + (j1 - j0 + 1) + " " + dx + " " + dy + "\n");
			for(int j = 0; j < data.length; j++) {
				for(int i = 0; i < data[0].length; i++) {
					for(int k = 0; k < data[0][0].length; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// insert
		// inserts data from one data file to another
		else if(args[0].toLowerCase().equals("i")) {
			s++;
			String from = args[s++];	// file name for data to be inserted
			String to = args[s++];		// file name for data where to be inserted
			String output = args[s++];

			// insertion start coordinates
			int i0 = 0;
			int j0 = 0;
			int k0 = 0;
			if(s < args.length) {
				i0 = s2i(args[s++]);
				j0 = s2i(args[s++]);
				k0 = s2i(args[s++]);
			}
			if(i0 < 0 || j0 < 0 || k0 < 0) {
				System.out.println("Error: Invalid insertion.");
				exit(1);
			}

			int Nx0 = Nx(from);		// source dimensions
			int Ny0 = Ny(from);
			int cols0 = cols(from);
			int Nx = Nx(to);		// destination dimensions
			int Ny = Ny(to);
			double dx = dx(to);
			double dy = dy(to);
			int cols = cols(to);
			int Nx_ = Math.max(Nx, i0 + Nx0);	// new dimensions
			int Ny_ = Math.max(Ny, j0 + Ny0);
			int cols_ = Math.max(cols, k0 + cols0);
			double[][][] data = new double[Ny_][Nx_][cols_];
			String[] words;

			// destination loaded first
			BufferedReader reader = read(to);
			reader.readLine();
			for(int j = 0; j < Ny; j++) {
				for(int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for(int k = 0; k < cols; k++) {
						data[j][i][k] = s2d(words[k]);
					}
				}
			}
			reader.close();

			// source loaded next (may overwrite destination data)
			reader = read(from);
			reader.readLine();
			for(int j = 0; j < Ny0; j++) {
				for(int i = 0; i < Nx0; i++) {
					words = reader.readLine().trim().split("\\s+");
					for(int k = 0; k < cols0; k++) {
						data[j0 + j][i0 + i][k0 + k] = s2d(words[k]);
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(Nx_ + " " + Ny_ + " " + dx + " " + dy + "\n");
			for(int j = 0; j < Ny_; j++) {
				for(int i = 0; i < Nx_; i++) {
					for(int k = 0; k < cols_; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// rotate
		// rotates data by 90, 180 or 270 degrees
		else if(args[0].toLowerCase().equals("r")) {
			s++;
			String input = args[s++];
			String output = args[s++];
			int phi = 90;				// rotation (default 90 degrees)
			if(s < args.length) {
				phi = s2i(args[s++]);	// rotation specified
			}

			int Nx = Nx(input);
			int Ny = Ny(input);
			double dx = dx(input);
			double dy = dy(input);

			String[][] data = new String[Ny][Nx];
			BufferedReader reader = read(input);
			reader.readLine();
			for(int j = 0; j < Ny; j++) {
				for(int i = 0; i < Nx; i++) {
					data[j][i] = reader.readLine().trim();
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			if(phi == 90 || phi == -270) {
				writer.write(Ny + " " + Nx + " " + dy + " " + dx + "\n");
				for(int i = 0; i < Nx; i++) {
					for (int j = Ny - 1; j >= 0; j--) {
						writer.write(data[j][i] + "\n");
					}
				}
			}
			else if(phi == 180 || phi == -180) {
				writer.write(Nx + " " + Ny + " " + dx + " " + dy + "\n");
				for(int j = Ny - 1; j >= 0; j--) {
					for(int i = Nx - 1; i >= 0; i--) {
						writer.write(data[j][i] + "\n");
					}
				}
			}
			else if(phi == 270 || phi == -90) {
				writer.write(Ny + " " + Nx + " " + dy + " " + dx + "\n");
				for(int i = Nx - 1; i >= 0; i--) {
					for(int j = 0; j < Ny; j++) {
						writer.write(data[j][i] + "\n");
					}
				}
			}
			else {
				System.out.println("Error: Invalid rotation.");
				exit(1);
			}
			writer.close();

		}

		// flip
		// flips data horizontally or vertically
		else if(args[0].toLowerCase().equals("f")) {
			s++;
			String input = args[s++];
			String output = args[s++];
			String direction = "x";						// flip direction
			if(s < args.length) direction = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			String[][] data = new String[Ny][Nx];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for(int j = 0; j < Ny; j++) {
				for(int i = 0; i < Nx; i++) {
					data[j][i] = reader.readLine().trim();
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			if(direction.equals("x") || direction.equals("h")) {	// x coordinates flipped
				for(int j = 0; j < Ny; j++) {
					for(int i = Nx - 1; i >= 0; i--) {
						writer.write(data[j][i] + "\n");
					}
				}
			}
			if(direction.equals("y") || direction.equals("v")) {	// y
				for(int j = Ny - 1; j >= 0; j--) {
					for(int i = 0; i < Nx; i++) {
						writer.write(data[j][i] + "\n");
					}
				}
			}
			writer.close();

		}

		// add
		// adds a constant value to data or adds to data together
		else if(args[0].toLowerCase().equals("add")) {
			s++;
			String input = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			double d;
			try {	// next argument a number? --> constant added

				d = Double.parseDouble(args[s]);
				s++;
				String output = args[s++];

				BufferedReader reader = read(input);
				String header = reader.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for (int k = 0; k < words.length; k++) {
							data[j][i][k] = s2d(words[k]) + d;
						}
					}
				}
				reader.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			} catch (NumberFormatException e) {	// not a number? --> another data set added

				String input_ = args[s++];
				String output = args[s++];

				int Nx_ = Nx(input_);
				int Ny_ = Ny(input_);
				int cols_ = cols(input_);

				if (Nx != Nx_ || Ny != Ny_ || cols != cols_) {	// dimensions match?
					System.out.println("Error: Incompatible dimensions.");
					exit(1);
				}

				String[] words_;
				BufferedReader reader = read(input);
				BufferedReader reader_ = read(input_);
				String header = reader.readLine();
				reader_.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						words_ = reader_.readLine().trim().split("\\s+");
						for (int k = 0; k < cols; k++) {
							data[j][i][k] = s2d(words[k]) + s2d(words_[k]);
						}
					}
				}
				reader.close();
				reader_.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			}
		}

		// subtract
		// similar to "add"
		else if(args[0].toLowerCase().equals("sub")) {
			s++;
			String input = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			double d;
			try {

				d = Double.parseDouble(args[s]);
				s++;
				String output = args[s++];

				BufferedReader reader = read(input);
				String header = reader.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for (int k = 0; k < words.length; k++) {
							data[j][i][k] = s2d(words[k]) - d;
						}
					}
				}
				reader.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			} catch (NumberFormatException e) {

				String input_ = args[s++];
				String output = args[s++];

				int Nx_ = Nx(input_);
				int Ny_ = Ny(input_);
				int cols_ = cols(input_);

				if (Nx != Nx_ || Ny != Ny_ || cols != cols_) {
					System.out.println("Error: Incompatible dimensions.");
					exit(1);
				}

				String[] words_;
				BufferedReader reader = read(input);
				BufferedReader reader_ = read(input_);
				String header = reader.readLine();
				reader_.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						words_ = reader_.readLine().trim().split("\\s+");
						for (int k = 0; k < cols; k++) {
							data[j][i][k] = s2d(words[k]) - s2d(words_[k]);
						}
					}
				}
				reader.close();
				reader_.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			}
		}

		// multiply
		// similar to "add"
		else if(args[0].toLowerCase().equals("mul")) {
			s++;
			String input = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			double d;
			try {

				d = Double.parseDouble(args[s]);
				s++;
				String output = args[s++];

				BufferedReader reader = read(input);
				String header = reader.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for (int k = 0; k < words.length; k++) {
							data[j][i][k] = s2d(words[k])*d;
						}
					}
				}
				reader.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			} catch (NumberFormatException e) {

				String input_ = args[s++];
				String output = args[s++];

				int Nx_ = Nx(input_);
				int Ny_ = Ny(input_);
				int cols_ = cols(input_);

				if (Nx != Nx_ || Ny != Ny_ || cols != cols_) {
					System.out.println("Error: Incompatible dimensions.");
					exit(1);
				}

				String[] words_;
				BufferedReader reader = read(input);
				BufferedReader reader_ = read(input_);
				String header = reader.readLine();
				reader_.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						words_ = reader_.readLine().trim().split("\\s+");
						for (int k = 0; k < cols; k++) {
							data[j][i][k] = s2d(words[k])*s2d(words_[k]);
						}
					}
				}
				reader.close();
				reader_.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			}
		}

		// divide
		// similar to "add"
		else if(args[0].toLowerCase().equals("div")) {
			s++;
			String input = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			double d;
			try {

				d = Double.parseDouble(args[s]);
				s++;
				String output = args[s++];

				BufferedReader reader = read(input);
				String header = reader.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for (int k = 0; k < words.length; k++) {
							data[j][i][k] = s2d(words[k])/d;
						}
					}
				}
				reader.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			} catch (NumberFormatException e) {

				String input_ = args[s++];
				String output = args[s++];

				int Nx_ = Nx(input_);
				int Ny_ = Ny(input_);
				int cols_ = cols(input_);

				if (Nx != Nx_ || Ny != Ny_ || cols != cols_) {
					System.out.println("Error: Incompatible dimensions.");
					exit(1);
				}

				String[] words_;
				BufferedReader reader = read(input);
				BufferedReader reader_ = read(input_);
				String header = reader.readLine();
				reader_.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						words_ = reader_.readLine().trim().split("\\s+");
						for (int k = 0; k < cols; k++) {
							data[j][i][k] = s2d(words[k])/s2d(words_[k]);
						}
					}
				}
				reader.close();
				reader_.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			}
		}

		// minimum
		// similar to "add"
		else if(args[0].toLowerCase().equals("min")) {
			s++;
			String input = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			double d;
			try {

				d = Double.parseDouble(args[s]);
				s++;
				String output = args[s++];

				BufferedReader reader = read(input);
				String header = reader.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for (int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.min(s2d(words[k]), d);
						}
					}
				}
				reader.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			} catch (NumberFormatException e) {

				String input_ = args[s++];
				String output = args[s++];

				int Nx_ = Nx(input_);
				int Ny_ = Ny(input_);
				int cols_ = cols(input_);

				if (Nx != Nx_ || Ny != Ny_ || cols != cols_) {
					System.out.println("Error: Incompatible dimensions.");
					exit(1);
				}

				String[] words_;
				BufferedReader reader = read(input);
				BufferedReader reader_ = read(input_);
				String header = reader.readLine();
				reader_.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						words_ = reader_.readLine().trim().split("\\s+");
						for (int k = 0; k < cols; k++) {
							data[j][i][k] = Math.min(s2d(words[k]), s2d(words_[k]));
						}
					}
				}
				reader.close();
				reader_.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			}
		}

		// maximum
		// similar to "add"
		else if(args[0].toLowerCase().equals("max")) {
			s++;
			String input = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			double d;
			try {

				d = Double.parseDouble(args[s]);
				s++;
				String output = args[s++];

				BufferedReader reader = read(input);
				String header = reader.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for (int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.max(s2d(words[k]), d);
						}
					}
				}
				reader.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			} catch (NumberFormatException e) {

				String input_ = args[s++];
				String output = args[s++];

				int Nx_ = Nx(input_);
				int Ny_ = Ny(input_);
				int cols_ = cols(input_);

				if (Nx != Nx_ || Ny != Ny_ || cols != cols_) {
					System.out.println("Error: Incompatible dimensions.");
					exit(1);
				}

				String[] words_;
				BufferedReader reader = read(input);
				BufferedReader reader_ = read(input_);
				String header = reader.readLine();
				reader_.readLine();
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						words_ = reader_.readLine().trim().split("\\s+");
						for (int k = 0; k < cols; k++) {
							data[j][i][k] = Math.max(s2d(words[k]), s2d(words_[k]));
						}
					}
				}
				reader.close();
				reader_.close();

				BufferedWriter writer = write(output);
				writer.write(header + "\n");
				for (int j = 0; j < Ny; j++) {
					for (int i = 0; i < Nx; i++) {
						for (int k = 0; k < cols; k++) {
							writer.write(data[j][i][k] + " ");
						}
						writer.write("\n");
					}
				}
				writer.close();

			}
		}

		// modulo
		// returns the modulo of data with given constant
		else if(args[0].toLowerCase().equals("mod")) {
			s++;
			String input = args[s++];
			double d = s2d(args[s++]);
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = s2d(words[k])%d;
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// power
		// similar to "mod"
		else if(args[0].toLowerCase().equals("pow")) {
			s++;
			String input = args[s++];
			double d = s2d(args[s++]);
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			if(d == 2.0) {
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for(int k = 0; k < words.length; k++) {
							d = s2d(words[k]);
							data[j][i][k] = d*d;
						}
					}
				}
			}
			else if(d == 0.5) {
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for(int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.sqrt(s2d(words[k]));
						}
					}
				}
			}
			else {
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for(int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.pow(s2d(words[k]), d);
						}
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// exponential
		// base can be specified, otherwise similar to "mod"
		else if(args[0].toLowerCase().equals("exp")) {
			s++;
			String input = args[s++];
			double base = Math.E;
			if(s < args.length - 1) {
				base = s2d(args[s++]);
			}
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			if(base == Math.E) {	// e^data?
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for(int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.exp(s2d(words[k]));
						}
					}
				}
			}
			else {					// base^data
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for(int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.pow(base, s2d(words[k]));
						}
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// logarithm
		// similar to "exp"
		else if(args[0].toLowerCase().equals("log")) {
			s++;
			String input = args[s++];
			double base = Math.E;
			if(s < args.length - 1) {
				base = s2d(args[s++]);
			}
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			if(base == Math.E) {
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for(int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.log(s2d(words[k]));
						}
					}
				}
			}
			else {
				double _lb = 1.0/Math.log(base);
				for(int j = 0; j < Ny; j++) {
					for(int i = 0; i < Nx; i++) {
						words = reader.readLine().trim().split("\\s+");
						for(int k = 0; k < words.length; k++) {
							data[j][i][k] = Math.log(s2d(words[k]))*_lb;
						}
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// sine
		// returns the sine of the data
		else if(args[0].toLowerCase().equals("sin")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.sin(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// cosine
		// similar to "sin"
		else if(args[0].toLowerCase().equals("cos")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.cos(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// tangent
		// ...
		else if(args[0].toLowerCase().equals("tan")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.tan(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// arcsine
		else if(args[0].toLowerCase().equals("asin")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.asin(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// arccosine
		else if(args[0].toLowerCase().equals("acos")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.acos(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// arctangent
		else if(args[0].toLowerCase().equals("atan")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.atan(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// hyperbolic sine
		else if(args[0].toLowerCase().equals("sinh")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.sinh(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// hyperbolic cosine
		else if(args[0].toLowerCase().equals("cosh")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.cosh(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// hyperbolic tangent
		else if(args[0].toLowerCase().equals("tanh")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						data[j][i][k] = Math.tanh(s2d(words[k]));
					}
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					for (int k = 0; k < cols; k++) {
						writer.write(data[j][i][k] + " ");
					}
					writer.write("\n");
				}
			}
			writer.close();

		}

		// complex conjugate
		// returns the complex conjugate of the data
		else if(args[0].toLowerCase().equals("conj")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			if(cols != 2) {		// not complex-valued data?
				System.out.println("Error: Invalid input. Complex-valued data expected.");
				exit(1);
			}

			String[] words;
			double[][][] data = new double[Ny][Nx][cols];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					data[j][i][0] = s2d(words[0]);
					data[j][i][1] = - s2d(words[1]);
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					writer.write(data[j][i][0] + " " + data[j][i][1] + "\n");
				}
			}
			writer.close();

		}

		// argument
		// similar to "conj"
		else if(args[0].toLowerCase().equals("arg")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);

			if(cols != 2) {
				System.out.println("Error: Invalid input. Complex-valued data expected.");
				exit(1);
			}

			String[] words;
			double[][] data = new double[Ny][Nx];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					data[j][i] = Math.atan2(s2d(words[1]), s2d(words[0]));
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					writer.write(data[j][i] + "\n");
				}
			}
			writer.close();

		}

		// magnitude
		// returns the Pythagorean magnitude of the data
		else if(args[0].toLowerCase().equals("abs")) {
			s++;
			String input = args[s++];
			String output = args[s++];

			int Nx = Nx(input);
			int Ny = Ny(input);
			int cols = cols(input);
			String[] words;
			double[][] data = new double[Ny][Nx];

			BufferedReader reader = read(input);
			String header = reader.readLine();
			double d;
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					words = reader.readLine().trim().split("\\s+");
					for (int k = 0; k < words.length; k++) {
						d = s2d(words[k]);
						data[j][i] += d*d;
					}
					data[j][i] = Math.sqrt(data[j][i]);
				}
			}
			reader.close();

			BufferedWriter writer = write(output);
			writer.write(header + "\n");
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {
					writer.write(data[j][i] + "\n");
				}
			}
			writer.close();

		}

		// garbage
		else {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

	}
}
