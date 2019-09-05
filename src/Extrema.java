// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package data;

import java.io.*;

import static java.lang.System.exit;

public class Extrema {

	public static void main(String[] args) throws IOException {

		// print instructions
		if(args.length == 0) {
			System.out.println("This is a tool for determining the minimum, maximum and average values of a real- or complex-valued field. The program expects the name of the input data file and optionally the name of the output text file:");
			System.out.println();
			System.out.println("java -jar extrema.jar image.png extrema.txt");
			exit(0);
		}

		if(args.length > 2) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

		String filename = args[0];
		double min = 1.0e300;
		double ave = 0.0;
		double max = -1.0e300;

		File file = new File(filename);
		if(!file.exists()) {
			System.out.println("Error: Input file \"" + filename + "\" not found.");
			exit(1);
		}

		BufferedReader reader = new BufferedReader(new FileReader(filename));
		String line = reader.readLine();
		String[] words = line.trim().split("\\s+");
		if(words.length != 4) {
			System.out.println("Error: Invalid input file header.");
			exit(1);
		}
		// system dimensions
		int Nx = Integer.parseInt(words[0]);
		int Ny = Integer.parseInt(words[1]);
		line = reader.readLine();
		words = line.trim().split("\\s+");
		double a = 0.0;
		double b = 0.0;
		double abs;
		// real-valued data
		if(words.length == 1) {
			try {
				a = Double.parseDouble(line);
			} catch (NumberFormatException e) {
				System.out.println("Error: Invalid input file. Decimal number expected instead of \"" + line + "\".");
				exit(1);
			}
			if(a < min) min = a;
			ave += a;
			if(a > max) max = a;
			while((line = reader.readLine()) != null) {
				try {
					a = Double.parseDouble(line);
				} catch (NumberFormatException e) {
					System.out.println("Error: Invalid input file. Decimal number expected instead of \"" + line + "\".");
					exit(1);
				}
				if(a < min) min = a;
				ave += a;
				if(a > max) max = a;
			}
			ave *= 1.0/Nx/Ny;
		}
		// complex-valued data
		else if(words.length == 2) {
			try {
				a = Double.parseDouble(words[0]);
			} catch (NumberFormatException e) {
				System.out.println("Error: Invalid input file. Decimal number expected instead of \"" + line + "\".");
				exit(1);
			}
			try {
				b = Double.parseDouble(words[1]);
			} catch (NumberFormatException e) {
				System.out.println("Error: Invalid input file. Decimal number expected instead of \"" + line + "\".");
				exit(1);
			}
			abs = Math.sqrt(a*a + b*b);
			if(abs < min) min = abs;
			ave += abs;
			if(abs > max) max = abs;
			while((line = reader.readLine()) != null) {
				words = line.trim().split("\\s+");
				try {
					a = Double.parseDouble(words[0]);
				} catch (NumberFormatException e) {
					System.out.println("Error: Invalid input file. Decimal number expected instead of \"" + line + "\".");
					exit(1);
				}
				try {
					b = Double.parseDouble(words[1]);
				} catch (NumberFormatException e) {
					System.out.println("Error: Invalid input file. Decimal number expected instead of \"" + line + "\".");
					exit(1);
				}
				abs = Math.sqrt(a*a + b*b);
				if(abs < min) min = abs;
				ave += abs;
				if(abs > max) max = abs;
			}
			ave *= 1.0/Nx/Ny;
		}
		else {
			System.out.println("Error: Invalid input file.");
			exit(1);
		}
		reader.close();

		System.out.println(min + " " + ave + " " + max);

		// write optional output file
		if(args.length == 2) {
			BufferedWriter writer = new BufferedWriter(new FileWriter(args[1]));
			writer.write(min + " " + ave + " " + max);
			writer.close();
		}

		exit(0);
	}
}
