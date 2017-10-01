package opt;

import pot.Benzene;
import pot.LJ;
import pot.Potential;
import pot.TIP4P;
import rotation.Rotation;
import rotation.SXNA;
import gop.Hop;
import line.LineSearch;
import line.WeakWolfe;
import line.Wolfe;

@SuppressWarnings("unused")
public class ConjGrad extends Optimiser {

	public double[] optimise(double[] x, LineSearch lineSearch, Potential pot, Rotation rotation) {
		
		double[] g = pot.grad(x, rotation);
		double[] gOld = new double[g.length];
		double beta = 0.0;
		double g1 = 0.0;
		double g2 = 0.0;
		
		double[] d = new double[g.length]; 
		double d1 = 0.0;
		
		double rms = 1.0d;
		double counter;
		
		double a = 0.0;
		double e = 1E-4;
		
		int iter = 0;
		
		double energy = pot.energy(x, rotation);
		double oldEnergy = 0;
		
		while(rms > e) {
			
			counter = 0;
			
			if(iter == 0) {
				
				for(int i = 0; i < d.length; i++) {
					d[i] = -g[i];
				}
				
				a = lineSearch.lS(1E-01, 0.1, g, x, d, pot, rotation);
				//System.out.println(a);
				
				for(int i = 0; i < x.length; i++) {
					x[i] = x[i] + a*d[i];
					counter += g[i]*g[i];
				}
			} else {
				for(int i = 0; i < x.length; i++) {
					g1 += g[i]*(g[i] - gOld[i]);
					g2 += gOld[i]*gOld[i];
				}
				
				beta = g1/g2;
				
				for(int i = 0; i < x.length; i++) {
					d[i] = -g[i] + beta*d[i];
				}
				
				a = lineSearch.lS(1E-01, 0.4, g, x, d, pot, rotation);
				//System.out.println(a);
				
				for(int i = 0; i < x.length; i++) {
					x[i] = x[i] + a*d[i];
					counter += g[i]*g[i];
				}
			}
			
			for(int i = 0; i < g.length; i++) {
				gOld[i] = g[i];
			}
			
			//System.out.println("ENERGY DIFF: " + ((energy-oldEnergy)==0.0));

			g = pot.grad(x, rotation);
			g1 = 0.0;
			g2 = 0.0;
			
			rms = (Math.sqrt(counter/g.length));
			iter++;
			
			oldEnergy = energy;
			energy = pot.energy(x, rotation);
			if((energy-oldEnergy) == 0.0) {
				rms = -1;
			}
			
			if(iter%50000 == 0) {
				System.out.println("FAILED");
				Hop.fail = true;
				break;
			}
			//System.out.println("ENERGY: " + energy);
		}
		System.out.println("ENERGY: " + energy);
		System.out.println("i: " + iter);
		Hop.iter(iter);
		return x;
	}
	
	public static void main(String[] args) {
		 
		 Potential pot = new TIP4P();
		 Benzene b = new Benzene();
		 
		 LineSearch w = new Wolfe();
		 ConjGrad cg = new ConjGrad();
		 Rotation sxna = new SXNA();
		 
		 double[] x = {-0.0804505953, 1.3318831021, 0.3396285110, 
				  	    0.0804505953,-3.3318831021,-0.3396285110,
				  	    0.8322660626,-0.7372539898, 0.1893163762,
				  	    1.8960194774, 1.3595002754, 0.1897039342};
		 
		 x = cg.optimise(x, w, b, sxna);
	}
}
