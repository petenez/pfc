// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package data;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

import static java.lang.System.exit;

public class IMG2PFC {

	public static void main(String[] args) throws IOException {
		if(args.length == 0) {
			System.out.println("This is a tool for converting image files into density fields. The average value and amplitude can be set using pfc-init.jar or manipulator.jar. An example of valid syntax:");
			System.out.println();
			System.out.println("java -jar img2n.jar image.png density.n");
			return;
		}

		if(args.length != 2) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}

		String input = args[0];
		String output = args[1];

		File file = new File(input);
		if(!file.exists()) {
			System.out.println("Error: Input file not found.");
			exit(1);
		}

		BufferedImage image = ImageIO.read(file);
		int W = image.getWidth();
		int H = image.getHeight();
		Writer writer = new BufferedWriter(new FileWriter(output));
		writer.write(W + " " + H + " 1.0 1.0\n");
		Color c;
		String out;
		for(int j = 0; j < H; j++) {
			for(int i = 0; i < W; i++) {
				c = new Color(image.getRGB(i, H - j - 1));
				out = "" + Color.RGBtoHSB(c.getRed(), c.getGreen(), c.getBlue(), null)[2];
				writer.write(out + "\n");
			}
		}
		writer.close();
	}
}