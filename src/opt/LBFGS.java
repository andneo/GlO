package opt;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import gop.Hop;
import line.LineSearch;
import line.WeakWolfe;
import line.Wolfe;
import pot.Benzene;
import pot.LJ;
import pot.Potential;
import pot.TIP4P;
import rotation.Rotation;
import rotation.SXNA;

@SuppressWarnings("unused")
public class LBFGS extends Optimiser{
	
	private BufferedWriter writer = null;

	public double[] optimise(double[] x, LineSearch lineSearch, Potential potential, Rotation rotation) {
		
		FilmWriter fw = new FilmWriter();
		
		double[] g = potential.grad(x, rotation);
		double[] gOld = new double[g.length]; 
		double[] d = new double[g.length]; 
		double gDot = 0.0;
		
		for(int i = 0; i < x.length; i++) {
			d[i] = -g[i];
		}
		
		double rms = 1.0d;
		
		double a = 0.0;
		double e = 1E-4;
		int iter = 0;
		
		int M = 4;
		int k = 0;
		int c = 0;
		int bound = 0;
		
		double[][] s = new double[M][x.length];
		double[][] y = new double[M][g.length];
		
		double p[] = new double[M];
		
		double alpha[] = new double[M];
		double beta = 0;
		double[] q = new double[g.length];
		double[] r = new double[d.length];
		double rk = 1.0;
		
		double ys = 1.0d;
		double yy = 1.0d;
		
		while(rms > e) {
			
			a = lineSearch.lS(1E-01, 0.4, g, x, d, potential, rotation);
			
			//System.out.println(a);
			for(int j = 0; j < x.length; j++) {
				x[j] += a*d[j];
				s[k][j] = a*d[j];
				gOld[j] = g[j];
				gDot += g[j]*g[j];
			}
			
			g = potential.grad(x, rotation);
			
			for(int j = 0; j < x.length; j++) {
				y[k][j] = g[j] - gOld[j];
				ys += y[k][j]*s[k][j];
				yy += y[k][j]*y[k][j];
				q[j] = g[j];
			}

			if(ys == 0.0) {
				ys = 1.0;
			}

			if(yy == 0.0) {
				yy = 1.0;
			}
			
			rk = ys/yy;
			p[k] = 1.0/ys;
			
			k++;
			
			if(k == M) {
				k = 0;
			}
				
			c = k;
			
			if(iter < M) {
				bound = iter;
			} else {
				bound = M;
			}
				
			for(int i = 0; i < bound; i++) {
				c -= 1;
				if(c == -1) {
					c = M - 1;
				}
				
				alpha[c] = 0.0;
				
				for(int j = 0; j < s[k].length; j++) {
					alpha[c] += p[c]*s[c][j]*q[j];
				}
				
				for(int j = 0; j < q.length; j++) {
					q[j] = q[j] - (alpha[c]*y[c][j]);
				}
			}
				
			for(int i = 0; i < r.length; i++) {
				r[i] = rk*q[i];
			}
				
			for(int i = 0; i < bound; i++) {
				
				beta = 0.0;
				
				for(int j = 0; j < s[k].length; j++) {
					beta += p[c]*y[c][j]*r[j];
				}
				
				for(int j = 0; j < r.length; j++) {
					r[j] = r[j] + s[c][j]*(alpha[c] - beta);
				}
				c += 1;
				if(c == M) {
					c = 0;
				}
			}
			
			for(int j = 0; j < d.length; j++) {
				d[j] = -r[j];
			}
			
			ys = 0.0d;
			yy = 0.0d;

			rms = (Math.sqrt(gDot/g.length));
			gDot = 0.0;
			iter++;
			
			//if(potential.toString().substring(4, 6).equals("LJ")) {
			//	fw.coordsFilm(x.length/3, x, "LJ", rotation, iter);
			//} else {
			//	fw.coordsFilm((x.length/2)/3, x, "TIP4P", rotation, iter);
			//}
			
			if(iter > 1000 || Double.isInfinite(potential.energy(x, rotation)) || Math.abs(potential.energy(x, rotation)) > 1E8) {
				System.out.println("FAILED");
				Hop.fail = true;
				rms = -1;
			}
			
			//System.out.println("ENERGY: " + potential.energy(x, rotation));

		}
		//System.out.println(potential.toString().substring(4, 7));
		//System.out.println("ENERGY: " + potential.energy(x, rotation));
		//System.out.println("i: " + iter);
		//System.out.println("\n");
		Hop.iter(iter);
		return x;
	}
	
	public static void main(String[] args) {
		 
		 LJ lj = new LJ();
		 TIP4P tip = new TIP4P();
		 Benzene b = new Benzene();
		 
		 LineSearch w = new Wolfe();
		 LBFGS lbfgs = new LBFGS();
		 Rotation sxna = new SXNA();
		 
		 double[] x = {-1.0804505953, 1.3318831021, 0.3396285110, 
				  	    1.0804505953,-3.3318831021,-0.3396285110,
				  	    0.8322660626,-0.7372539898, 0.1893163762,
				  	    1.8960194774, 1.3595002754, 0.1897039342};
		 
		 //double[] x = {1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,0,1, 0,1,1};
		 
		 x = lbfgs.optimise(x, w, b, sxna);
	}
}