package line;

import pot.Potential;
import rotation.Rotation;

/**
 * 
 * @author andreasneophytou
 * Projected search methods
 */

public class BackTrack extends LineSearch {
	
	public double lS(double c1, double t, double[] g, double[] x, double[] d, Potential pot, Rotation rotation) {
		double a = 1.0;
		
		c1 = 0.5;
		t = 0.5;		
		
		double f0 = pot.energy(x, rotation);
		double[] x1 = new double[x.length];
		double f_p0 = 0.0;
		double f_i = 0.0;
		
		for(int i = 0; i < x1.length; i++) {
			x1[i] = x[i] + a*d[i];
			f_p0 += g[i]*d[i];
		}
		
		f_i = pot.energy(x1, rotation);
		
		while(f_i > f0 + a*c1*f_p0) {
			a *= t;
			
			for(int i = 0; i < x1.length; i++) {
				x1[i] = x[i] + a*d[i];
			}
			
			f_i = pot.energy(x1, rotation);
		}
		
		//System.out.println("No. of Linesearch steps: " + it);
		
		if(f_i > f0) {
			//System.out.println("\nERROR, a: " + a);
			while(f_i > f0) {
				a = a*(1E-15);
				for(int i = 0; i < x1.length; i++) {
					x1[i] = x[i] + a*d[i];
				}
				f_i = pot.energy(x1, rotation);
			}
			
			//System.out.println("**********************************************************************************************");
			//System.out.println("Old Energy: " + pot.energy(x) + ", New Energy: " + pot.energy(x1) + ", a: " + a);
			//System.out.println("**********************************************************************************************");
		}
		
		return a;
	}

}
