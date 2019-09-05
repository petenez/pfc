// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package visualization;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static java.lang.System.exit;

public class AltPlotter {
	 
	public static void main(String[] args) throws IOException {

		if(args.length == 0) {
			System.out.println("This is an alternative version of plotter.jar. Real-valued fields are plotted with an inverted color scale (minimum white, maximum black). For complex-valued fields, the magnitude is mapped to saturation and the phase to hue. See plotter.jar for further instructions.");
			exit(0);
		}

		String input;
		String output;
		int Lx;	// width of the system in pixels
		int Ly;	// height
		// parameters for optional looping over time steps
		int start;
		int incr;
		int end;
		if(args.length != 2 && args.length != 5) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}
		if(args.length == 5) {
			start = Integer.parseInt(args[2]);
			incr = Integer.parseInt(args[3]);
			end = Integer.parseInt(args[4]);
			if(start > end || incr < 1) {
				System.out.println("Error: Invalid batch mode parameters.");
				exit(1);
			}
		} else {
			start = 0;
			incr = 1;
			end = 0;
		}
		
		boolean real;	// real-valued data?
		
		// loop over different time steps
		for(int m = start; m <= end; m += incr) {
			input = args[0];
			output = args[1];
			if(args.length == 5) {	// loop over different time steps
				String[] pieces = input.trim().split("#");	// split 'input' at "#"
				input = "";	// input emptied
				if(pieces.length >= 1) input += pieces[0];	// append piece before "#"
				input += m;	// append time step
				if(pieces.length >= 2) input += pieces[1];	// append piece after "#" 
				pieces = output.trim().split("#");	// repeat for 'output'
				output = "";
				if(pieces.length >= 1) output += pieces[0];
				output += m;
				if(pieces.length >= 2) output += pieces[1];
			}
			File f = new File(input);
			if(!f.exists()) {	// check if input file is found
				if(start < end) {
					System.out.println("Warning: Input file " + input + " not found. Continuing iteration.");
					continue;
				}
				else {
					System.out.println("Error: Input file " + input + " not found.");
					exit(1);
				}
			}
			
			// checks if the data is real- or complex-valued
			BufferedReader reader = new BufferedReader(new FileReader(input));
			reader.readLine();
			int parts = reader.readLine().trim().split("\\s+").length;	// splits the line at " "
			if(parts == 1) real = true;	// single column -> real
			else if(parts == 2) real = false;	// double -> complex
			else {
				real = false;
				System.out.println("Error: Invalid data.");
				reader.close();
				exit(1);
			}
			reader.close();
			reader = new BufferedReader(new FileReader(input));	// start again
			String[] temp = reader.readLine().trim().split("\\s+");
			Lx = Integer.parseInt(temp[0]);
			Ly = Integer.parseInt(temp[1]);
			double[][][] data = new double[Lx][Ly][2];	// data array

			// determines extrema
			double min = 1.0e300;	// initially huge just to be safe
			double max = -1.0e300;	// initially -huge ...
			if(real) {
				double a;
				for(int j = 0; j < Ly; j++) {
					for(int i = 0; i < Lx; i++) {
						a = Double.parseDouble(reader.readLine());
						data[i][j][0] = a;
						if(a < min) min = a;	// check if current minimum
						if(a > max) max = a;	// maximum
					}
				}
			} else {
				double a, b, A;
				min = 0.0;
				for(int j = 0; j < Ly; j++) {
					for(int i = 0; i < Lx; i++) {
						temp = reader.readLine().trim().split("\\s+");
						a = Double.parseDouble(temp[0]);
						data[i][j][0] = a;
						b = Double.parseDouble(temp[1]);
						data[i][j][1] = b;
						A = Math.sqrt(a*a+b*b);	// magnitude
						if(A > max) max = A;
					}
				}
			}
			reader.close();

			// sets pixels and writes output
			BufferedImage image = new BufferedImage(Lx, Ly, BufferedImage.TYPE_INT_RGB);
			double a, b, A;
			float arg;
			if(real) {	// real-valued data
				for(int i = 0; i < Lx; i++) {
					for(int j = 0; j < Ly; j++) {
						image.setRGB(i, Ly-j-1, Color.HSBtoRGB(0, 0, (float)((max - data[i][j][0])/(max - min))));
					}
				}
			} else {	// complex-valued data
				double div2pi = 0.5/Math.PI;
				for(int i = 0; i < Lx; i++) {
					for(int j = 0; j < Ly; j++) {
						a = data[i][j][0];
						b = data[i][j][1];
						A = Math.sqrt(a*a + b*b);
						arg = (float)(div2pi*Math.atan2(b, a) + 0.5);	// phase from [-pi, pi] to [0, 1]
						image.setRGB(i, Ly - j - 1, Color.HSBtoRGB(arg, (float)(A/(max - min)), 1));
					}
				}
			}
			ImageIO.write(image, "png", new File(output));
		}
	}
}
