package line;

import pot.Potential;
import rotation.Rotation;

/**
 * 
 * @author andreasneophytou
 * Linesearch using Armijo and Wolfe criteria. 
 */

public class Wolfe extends LineSearch {
	
	public double lS(double c1, double c2, double[] g, double[] x, double[] d, Potential pot, Rotation rotation) {
		
		double aOld = 0.0; // Set a_0 to be equal to zero
		double aHigh = 200.0; // Set aHigh to be arbitrarily large
		double a = 1.0; // Set a_1 to be equal to 1.0
	
		double f0 = pot.energy(x, rotation);
		
		double f_p0 = 0.0;
		for(int i = 0; i < x.length; i++) {
			f_p0 += g[i]*d[i]; 
		}
		
		double[] x_i = new double[x.length];
		double[] g_i = new double[x.length];
		
		double f_i = 0.0;
		double f_iOld = f0;
		double f_pi = 0.0;
		
		boolean done = false;
		
		int it = 0;
		
		while(done == false || a == aHigh) {
			
			for(int i = 0; i < x_i.length; i++) {
				x_i[i] = x[i] + a*d[i];
			}
			
			f_i = pot.energy(x_i, rotation);
			
			if(f_i > f0 + a*c1*(f_p0) || (it > 0 && f_i >= f_iOld)) {
				a = zoom(aOld,a,f0,f_p0,c1,c2,x,d,pot,rotation);
				done = true;
				break;
			} 
			
			f_pi = 0.0;
			g_i = pot.grad(x_i, rotation);
			for(int i = 0; i < x_i.length; i++) {
				f_pi += g_i[i]*d[i]; 
			}
			
			if(Math.abs(f_pi) <= -c2*f_p0) {
				done = true;
				break;
			} else if(f_pi >= 0) {
				a = zoom(a,aOld,f0,f_p0,c1,c2,x,d,pot,rotation);
				done = true;
				break;
			} 
			
			aOld = a;
			a *= 2.0;
			//a = (a+aHigh)/2.0;
			
			f_iOld = f_i;
			it++;
		}
		return a;
	}
	
	public static double zoom(double aLow, double aHigh, double f0, double f_p0, double c1, double c2, double[] x, double[] d, Potential pot, Rotation rotation) {
		
		double a = 0.0;
		
		double[] x_j = new double[x.length];
		double[] g_j = new double[x.length];
		double f_j = 0.0;
		double f_pj = 0.0;
		
		double f_Low = 0.0;
		double[] x_Low = new double[x.length];
		
		boolean done = false;
		
		int iter = 0;
		
		while(done == false) {
			
			a = (aLow+aHigh)/2.0;
		
			for(int i = 0; i < x_j.length; i++) {
				x_j[i] = x[i] + a*d[i];
				x_Low[i] = x[i] + aLow*d[i];
			}
		
			f_j = pot.energy(x_j, rotation);
			
			f_Low = pot.energy(x_Low, rotation);
		
			if(f_j > f0 + a*c1*f_p0|| f_j >= f_Low) {
				aHigh = a;
			} else {
				g_j = pot.grad(x_j, rotation);
				f_pj = 0.0;
				for(int i = 0; i < x_j.length; i++) {
					f_pj += g_j[i]*d[i]; 
				}
				if(Math.abs(f_pj) <= -c2*f_p0) {
					done = true;
					break;
				} else if(f_pj*(aHigh-aLow) >= 0) {
					aHigh = aLow;
				} 
			
				aLow = a;
			}
			
			if(iter == 500) {
				break;
			}
			
			iter++;
		}
		return a;		
	}
}
