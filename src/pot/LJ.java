package pot;

import rotation.Rotation;

public class LJ extends Potential {
	
	public int factorial(int n) {
		int x = 0;
		for(int i = 1; i < n; i++) {
			x += i;
		}
		return x;
	}
	
	public double[] dist(double[] x) {
		
		int natm = x.length/3;
		double[] dist = new double[factorial(natm)*3];
		int k = 0;
		
		for(int i = 0; i < natm - 1; i++) {
			for(int j = i + 1; j < natm; j++) {
				dist[k] = Math.abs(x[(i*3)] - x[(j*3)]);
				dist[k+1] = Math.abs(x[(i*3)+1] - x[(j*3)+1]);
				dist[k+2] = Math.abs(x[(i*3)+2] - x[(j*3)+2]);
				k += 3;
			}
		}
		
		return dist;
	}
	
	public double energy(double[] x, Rotation rotation) {
		
		int natm = x.length/3;
		double energy = 0;
		double[] Rij = new double[3];
		double rij = 0.0;
		
		for(int i = 0; i < natm - 1; i++) {
			for(int j = i + 1; j < natm; j++) {
				
				Rij[0] = x[(i*3)] - x[(j*3)];
				Rij[1] = x[(i*3)+1] - x[(j*3)+1];
				Rij[2] = x[(i*3)+2] - x[(j*3)+2];
				
				rij = Math.sqrt((Rij[0]*Rij[0])+(Rij[1]*Rij[1])+(Rij[2]*Rij[2]));
				energy += 4*(Math.pow(rij, -12) - Math.pow(rij, -6));
			}
		}
		
		return energy;
	}
	
	public double[] grad(double[] x, Rotation rotation) {
		
		int natm = x.length/3;
		double[] grad = new double[x.length];
		
		double[] Rij = new double[3];
		double rij = 0.0;
		double gE = 0.0;
		
		for(int i = 0; i < natm - 1; i++) {
			for(int j = i + 1; j < natm; j++) {
				
				Rij[0] = x[(i*3)] - x[(j*3)];
				Rij[1] = x[(i*3)+1] - x[(j*3)+1];
				Rij[2] = x[(i*3)+2] - x[(j*3)+2];	
				
				rij = Math.sqrt((Rij[0]*Rij[0])+(Rij[1]*Rij[1])+(Rij[2]*Rij[2]));
				gE = (Math.pow(rij, -8))*((2*Math.pow(rij, -6)) - 1);		
				
				grad[(i*3)] 	-= 24*(Rij[0]*gE);
				grad[(i*3)+1] 	-= 24*(Rij[1]*gE);
				grad[(i*3)+2] 	-= 24*(Rij[2]*gE);
				
				grad[(j*3)] 	+= 24*(Rij[0]*gE);
				grad[(j*3)+1] 	+= 24*(Rij[1]*gE);
				grad[(j*3)+2] 	+= 24*(Rij[2]*gE);
			}	
		}
		return grad;
	}
}
