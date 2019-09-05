// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package atomistic;

import java.util.ArrayList;

public class Atom {	// implements atoms and rings

	private String type;
	private int k;		// index of atom/ring
	private double[] r;	// its position
	// for an atom, its bros are its neighboring atoms
	// for a ring, its bros are its neighboring rings
	private ArrayList<Atom> bros = new ArrayList<Atom>();
	// for an atom, its hos are its neighboring rings
	// for a ring, its hos are its neighboring atoms
	private ArrayList<Atom> hos = new ArrayList<Atom>();

	public Atom(int k, double[] r) {
		this.k = k;
		this.r = r;
	}

	public String type() {
		return type;
	}

	public void type(String str) {
		type = str;
	}

	public int k() {
		return k;
	}

	public double[] r() {
		return r;
	}

	public void r(double[] r) {
		this.r = r;
	}

	public void add_bro(Atom p) {
		bros.add(p);
	}

	public Atom get_bro(int n) {
		return bros.get(n);
	}

	public boolean is_bro(Atom p) {
		if(bros.contains(p)) return true;
		return false;
	}

	public int count_bros() {
		return bros.size();
	}

	public void add_ho(Atom p) {
		hos.add(p);
	}

	public Atom get_ho(int n) {
		return hos.get(n);
	}

	public boolean is_ho(Atom p) {
		if(hos.contains(p)) return true;
		return false;
	}

	public int count_hos() {
		return hos.size();
	}
}