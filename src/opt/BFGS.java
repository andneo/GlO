package opt;

import gop.Hop;
import line.LineSearch;
import pot.Potential;
import rotation.Rotation;


public class BFGS extends Optimiser {

	public double[] optimise(double[] x, LineSearch lineSearch, Potential pot, Rotation rotation) {
		
		double[] g = pot.grad(x, rotation);
		double[] gOld = new double[g.length]; 
		double gDot = 0.0;
		
		double[] d = new double[g.length];
		
		double[][] C = new double[x.length][x.length];
		double[][] D = new double[x.length][x.length];
		
		for(int i = 0; i < x.length; i++) {
			for(int j = 0; j < x.length; j++) {
				C[i][j] = 0.0;
				D[i][j] = 0.0;
			}
		}
		
		double[][] H = new double[x.length][x.length];
		double[][] V = new double[x.length][x.length];
		
		for(int i = 0; i < x.length; i++) {
			for(int j = 0; j < x.length; j++) {
				if(i == j) {
					H[i][j] = 1.0;
					V[i][j] = 1.0;
				} else {
					H[i][j] = 0.0;
					V[i][j] = 0.0;
				}
			}
		}
		
		double rms = 1.0d;
		
		double a = 0.0;
		double e = 1E-4;
		
		int iter = 0;
		
		double[] s = new double[x.length];
		double[] y = new double[g.length];
		
		double p = 0.0;

		double ys = 0.0;
		
		while(rms > e) {
			
			for(int i = 0; i < H.length; i++) {
				d[i] = 0;
				for(int j = 0; j < H[0].length; j++) {
					d[i] -= H[i][j]*g[j];
				}
			}
			
			
			a = lineSearch.lS(1E-01, 0.4, g, x, d,pot, rotation);
			
			for(int i = 0; i < x.length; i++) {
				x[i] += a*d[i];
				s[i] = a*d[i];
				gOld[i] = g[i];
				gDot += g[i]*g[i];
			}
			
			g = pot.grad(x, rotation);
			
			for(int i = 0; i < x.length; i++) {
				y[i] = g[i] - gOld[i];
				ys += y[i]*s[i];
			}
			
			p = 1.0/ys;
			
			for(int i = 0; i < s.length; i++) {
				for(int j = 0; j < s.length; j++) {
					if(i == j) {
						V[i][j] = 1.0 - (p*y[i]*s[j]);
					} else {
						V[i][j] = -(p*y[i]*s[j]);
					}
				}
			}
			
			for(int i = 0; i < H.length; i++) {
				for(int j = 0; j < H[0].length; j++) {
					for(int k = 0; k < H[0].length; k++) {
						C[i][j] += H[i][k]*V[k][j];
					}
				}
			}
			
			for(int i = 0; i < H.length; i++) {
				for(int j = 0; j < H[0].length; j++) {
					for(int k = 0; k < H[0].length; k++) {
						D[i][j] += V[k][i]*C[k][j];
					}
				}
			}
			
			for(int i = 0; i < H.length; i++) {
				for(int j = 0; j < H[0].length; j++) {
					H[i][j] = D[i][j]+(p*s[i]*s[j]);
				}
			}
			
			ys = 0.0;
			
			for(int i = 0; i < C.length; i++) {
				for(int j = 0; j < C[0].length; j++) {
					C[i][j] = 0.0;
					D[i][j] = 0.0;
				}
			}
			
			rms = (Math.sqrt(gDot/g.length));
			
			gDot = 0.0;
			iter++;
			//System.out.println("ENERGY: " + pot.energy(x));
			if(iter > 50000) {
				rms = -1;
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
