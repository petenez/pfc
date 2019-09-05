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

public class Orienter {

	public static void main(String[] args) throws IOException {

		if(args.length == 0) {
			System.out.println("This is a tool for plotting an \"oriented\" PFC density field. This means that a PFC density field is colored based on the corresponding orientation field (the fields tool can be used to obtain the latter). In the final image, hue is given by the phase of the orientation field, saturation is given by the magnitude of the orientation field and brightness is given by the magnitude of the density field.");
			System.out.println();
			System.out.println("An example of valid syntax:");
			System.out.println();
			System.out.println("java -jar orienter.jar density.n orientation.phi result.png");
			exit(0);
		}

		if(args.length != 3) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

		String input1 = args[0];
		String input2 = args[1];
		String output = args[2];

		if(!(new File(input1)).exists() || !(new File(input2)).exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}

		BufferedReader reader1 = new BufferedReader(new FileReader(input1));
		BufferedReader reader2 = new BufferedReader(new FileReader(input2));

		String[] line = reader1.readLine().trim().split("\\s+");
		int Nx = Integer.parseInt(line[0]);
		int Ny = Integer.parseInt(line[1]);
		line = reader2.readLine().trim().split("\\s+");
		if(Integer.parseInt(line[0]) != Nx || Integer.parseInt(line[1]) != Ny) {
			System.out.println("Error: Mismatching dimensions.");
			exit(1);
		}

		BufferedReader reader = new BufferedReader(new FileReader(input1));
		reader.readLine();
		line = reader.readLine().trim().split("\\s+");
		if(line.length != 1) {
			System.out.println("Error: Invalid data.");
			reader.close();
			exit(1);
		}
		reader.close();
		reader = new BufferedReader(new FileReader(input2));
		reader.readLine();
		line = reader.readLine().trim().split("\\s+");
		if(line.length != 2) {
			System.out.println("Error: Invalid data.");
			reader.close();
			exit(1);
		}
		reader.close();

		double min1 = 1.0e300;
		double max1 = -1.0e300;
		double max2 = -1.0e300;
		double re, im, abs;
		double[][] n = new double[Ny][Nx];
		double[][][] phi = new double[Ny][Nx][2];
		for(int j = 0; j < Ny; j++) {
			for(int i = 0; i < Nx; i++) {
				n[j][i] = Double.parseDouble(reader1.readLine());
				if(n[j][i] < min1) min1 = n[j][i];
				if(n[j][i] > max1) max1 = n[j][i];
				line = reader2.readLine().trim().split("\\s+");
				phi[j][i][0] = Double.parseDouble(line[0]);
				phi[j][i][1] = Double.parseDouble(line[1]);
				re = phi[j][i][0];
				im = phi[j][i][1];
				abs = Math.sqrt(re*re + im*im);
				if(abs > max2) max2 = abs;
			}
		}

		BufferedImage image = new BufferedImage(Nx, Ny, BufferedImage.TYPE_INT_RGB);

		float h, s, b;
		double _2pi = 0.5/Math.PI;
		for(int j = 0; j < Ny; j++) {
			for(int i = 0; i < Nx; i++) {
				re = phi[j][i][0];
				im = phi[j][i][1];
				h = (float)(Math.atan2(im, re)*_2pi + 0.5);
				abs = Math.sqrt(re*re + im*im);
				s = (float)(abs/max2);
				b = (float)((n[j][i] - min1)/(max1 - min1));
				image.setRGB(i, Ny - j - 1, Color.HSBtoRGB(h, s, b));
			}
		}

		ImageIO.write(image, "png", new File(output));

		exit(0);
	}
}
