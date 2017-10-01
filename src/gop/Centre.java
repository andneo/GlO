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

}
