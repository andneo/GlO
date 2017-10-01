package gop;

import java.util.Comparator;

public class Minimum implements Comparator<Minimum> , Comparable<Minimum> {
	
	private double[] x;
	private double e;
	private int step;
	
	public double[] getX() {
		return x;
	}
	public void setX(double[] x) {
		this.x = x;
	}
	public double getE() {
		return e;
	}
	public void setE(double e) {
		this.e = e;
	}
	public int getStep() {
		return step;
	}
	public void setStep(int step) {
		this.step = step;
	}
	
	public Minimum(double[] x, double e, int step) {
		this.x = x;
		this.e = e;
		this.step = step;
	}
	
	public int compareTo(Minimum m) {
		return compare(this, m);
	}

	public int compare(Minimum m1, Minimum m2) {
		if(m1.e < m2.e) {
			return -1;
		} else {
			return 1;
		}
	}

}
