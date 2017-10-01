package main;

import gop.Hop;
import gui.View;

public class Main {
	
	public static void main(String[] args) {
		
		//new View();
		
		//String opt = "SD";
		//String opt = "CG";
		//String opt = "BFGS";
		String opt = "LBFGS";
		
		//String lS = "BackTrack";
		//String lS = "WeakWolfe";
		String lS = "Wolfe";
		
		String pot = "LJ";
		pot = "TIP";
		
		//String o = "AA";
		String o = "SXNA";
		
		Hop s = new Hop();
		
		//s.setParameters(int n, int nAtoms, double c, int seed, String LOpt, String lineSearch, String potential, String orientation, boolean random, double tStep, double rStep)
		s.gOptimise(51, 200, 10.0, 26, opt, lS, pot, o, false, 1.0, 1.5);
	}

}
