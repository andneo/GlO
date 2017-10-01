package line;

import pot.Potential;
import rotation.Rotation;

/**
 * 
 * @author andreasneophytou
 * Wolfe-Powell inexact line search. 
 * Implemented as described in Chapter 2 of:
 * Sun, Wenyu & Yuan, Ya-xiang. Optimization theory and methods. Nonlinear programming. (2006).
 */

public class WeakWolfe extends LineSearch {
	
	public double lS(double c1, double c2, double[] g, double[] x, double[] d, Potential pot, Rotation rotation) {
		double a = 1.0;
		double a1 = 0.0;
		double a2 = Double.MAX_VALUE;
	
		double f0 = pot.energy(x, rotation);
		
		double f_p0 = 0.0;
		for(int i = 0; i < x.length; i++) {
			f_p0 += g[i]*d[i]; 
		}
		
		double[] x_i = new double[x.length];
		double[] g_i = new double[x.length];
		
		double f_i = 0.0;
		double f_pi = 0.0;
		
		int iter = 0;
		
		while(true) {
			
			for(int i = 0; i < x_i.length; i++) {
				x_i[i] = x[i] + a*d[i];
			}
			
			f_i = pot.energy(x_i, rotation);
			g_i = pot.grad(x_i, rotation);
			
			for(int i = 0; i < x_i.length; i++) {
				f_pi += g_i[i]*d[i]; 
			}
			
			if(f_i > f0 + a*c1*(f_p0)) {
				a2 = a;
			} else if(f_pi < c2*f_p0) {
				a1 = a;
			} else {
				break;
			}
			
			if(a2 < Double.MAX_VALUE) {
				a = (a1+a2)/2.0;
			} else {
				a = 2.0*a;
			}	
			
			f_pi = 0.0;
			
			if(iter == 500) {
				break;
			}
			
			iter++;
		}
		
		return a;
	}
}
