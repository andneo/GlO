package opt;

import gop.Hop;
import line.LineSearch;
import pot.Potential;
import rotation.Rotation;

public class StDes extends Optimiser {

	public double[] optimise(double[] x, LineSearch lineSearch, Potential pot, Rotation rotation) {
		
		double[] g = pot.grad(x, rotation);
		double[] d = g;
		double rms = 1.0d;
		double counter;
		
		double a = 0.0;
		double e = 1E-4;
		
		int iter = 0;
		
		double energy = pot.energy(x, rotation);
		double oldEnergy = 0;
		
		while(rms > e) {
			
			counter = 0;
			
			g = (pot.grad(x, rotation));
			
			for(int i = 0; i < x.length; i++) {
				d[i] = -g[i];
			}
			
			a = lineSearch.lS(1E-01, 0.1, g, x, d,pot, rotation);
			
			for(int i = 0; i < x.length; i++) {
				x[i] = x[i] + a*d[i];
				counter += g[i]*g[i];
			}
			
			rms = (Math.sqrt(counter/g.length));
			iter++;
			
			oldEnergy = energy;
			energy = pot.energy(x, rotation);
			if((energy-oldEnergy) == 0.0) {
				rms = -1;
			}
			
			if(iter%500000 == 0) {
				System.out.println("FAILED");
				break;
			}
		}
		
		//System.out.println("ENERGY: " + pot.energy(x));
		//System.out.println("i: " + iter);
		Hop.iter(iter);
		return x;
	}
	
	public static void main(String[] args) {
	}
}