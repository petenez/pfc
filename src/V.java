// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package atomistic;

public class V {	// implements basic vector math

	public static double abs2(double[] a) {
		return a[0]*a[0]+a[1]*a[1];
	}

	public static double abs(double[] a) {
		return Math.sqrt(abs2(a));
	}

	public static double arg(double[] a) {
		return Math.atan2(a[1], a[0]);
	}

	public static double[] sum(double[] a, double[] b) {
		return new double[] {a[0] + b[0], a[1] + b[1]};
	}

	public static double[] sub(double[] a, double[] b) {
		return new double[] {a[0] - b[0], a[1] - b[1]};
	}

	public static double[] mul(double c, double[] a) {
		return new double[] {c*a[0], c*a[1]};
	}

	public static double[] div(double[] a, double c) {
		return new double[] {a[0]/c, a[1]/c};
	}

	// cross product
	public static double x(double[] a, double[] b) {
		return a[0]*b[1] - a[1]*b[0];
	}

	// cross product with out-of-plane unit vector
	public static double[] x(double[] a) {
		return new double[] {a[1], -a[0]};
	}

	// dot product
	public static double dot(double[] a, double[] b) {
		return a[0]*b[0] + a[1]*b[1];
	}

	// angle between a and b
	public static double angle(double[] a, double[] b) {
		return Math.acos(dot(a, b)/(abs(a)*abs(b)));
	}
}