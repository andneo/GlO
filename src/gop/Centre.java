package gop;

public class Centre {
	
	public double[] centre(double[] x, int nAtoms) {
		
		double cx = 0.0;
		double cy = 0.0;
		double cz = 0.0;
		int counter = 0;
		
		for(int j = 0; j < nAtoms; j++) {
			counter = (j+1)*3;
			cx += x[counter-3];
			cy += x[counter-2];
			cz += x[counter-1];
		}
		
		cx = cx/nAtoms;
		cy = cy/nAtoms;
		cz = cz/nAtoms;
		
		for(int j = 0; j < nAtoms; j++) {
			counter = (j+1)*3;
			x[counter-3] -= cx;
			x[counter-2] -= cy;
			x[counter-1] -= cz;
		}
		
		return x;
	}
	
	public boolean tooClose(double[] x) {
		boolean tClose = false;
		double[] Rij = new double[3];
		double rij = 0.0;
		
		for(int i = 0; i < x.length/2; i+=3) {
			for(int j = i+3; j < x.length/2; j+=3) {
				Rij[0] = x[i] - x[j];
				Rij[1] = x[i+1] - x[j+1];
				Rij[2] = x[i+2] - x[j+2];
				
				rij = Math.sqrt(Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2]); 
				
				if(rij < 2.3) {
					tClose = true;
					//System.out.println(rij + ", Too Close");
					System.out.println(Rij[0] + " " +Rij[1] + " " + Rij[2]);
					break;
				}
			}
			
			if(tClose) break;
		}
		
		return tClose;
	}

}
