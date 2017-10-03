package pot;

import rotation.Mat;
import rotation.Rotation;
import rotation.SXNA;

public class Benzene extends Potential {
	
	private double[][] rbsites = {{ 1.3970, 0.0000, 0.0000},
								  {-1.3970, 0.0000, 0.0000},
								  { 0.6985, 1.2098, 0.0000},
								  {-0.6985,-1.2098, 0.0000},
								  {-0.6985, 1.2098, 0.0000},
								  { 0.6985,-1.2098, 0.0000},
								  { 2.4240, 0.0000, 0.0000},
								  {-2.4240, 0.0000, 0.0000},
								  { 1.2120, 2.0992, 0.0000},
								  {-1.2120,-2.0992, 0.0000},
								  {-1.2120, 2.0992, 0.0000},
								  { 1.2120,-2.0992, 0.0000}};
	
	private double[] q = {-0.153, 0.153};
	
	private double[] A = {2414.0,	573.0,	 136.0};
	private double[] B = {367250.0,	65485.0, 11677.0};
	private double[] C = {3.60,		3.67,	 3.74};
	private double coulomb = 1389.354848;
	
	//private double[] mass = {12.0, 1.0};
	private double[] cm = new double[3];
	private double m = 0.0;
	
	private String orientation;
	
	public String getOrientation() {
		return orientation;
	}

	public void setOrientation(String orientation) {
		this.orientation = orientation;
	}

	private double[][] getRBSites() {
		
		for(int i = 0; i < 12; i++) {
			if(i < 6) {
				cm[0] += 12.0*rbsites[i][0]; 
				cm[1] += 12.0*rbsites[i][1]; 
				cm[2] += 12.0*rbsites[i][2];
				m 	  += 12.0;
			} else {
				cm[0] += rbsites[i][0]; 
				cm[1] += rbsites[i][1]; 
				cm[2] += rbsites[i][2];
				m 	  += 1.0;
			}
		}
		
		for(int i = 0; i < 3; i++) {
			cm[i] = cm[i]/m; 
		}
		
		for(int i = 0; i < 12; i++) {
			rbsites[i][0] -= cm[0];
			rbsites[i][1] -= cm[1];
			rbsites[i][2] -= cm[2];
		}
		
		return rbsites;	
	}
	
	public double[][] viewBZ(double[] x, Rotation rotation) {
		
		double[] p      	= new double[3];
		double[][] RM   	= new double[3][3];
		double[][] DRMX 	= new double[3][3];
		double[][] DRMY 	= new double[3][3];
		double[][] DRMZ 	= new double[3][3];
		double[] rb     	= new double[3];
		double[][] rbcoords = new double[((x.length/2)/3)*12][3];
		int offset 			= x.length/2;
		int k = 0;
		
		for(int i = 0; i < offset; i+=3) {
			p[0] = x[i+offset];
			p[1] = x[i+1+offset];
			p[2] = x[i+2+offset];
			
			Mat mat = new Mat();
			rotation.rotDrvt(p, RM, DRMX, DRMY, DRMZ, false);
			
			for(int j = 0; j < 12; j++) {
				rb = mat.vdot(RM,  getRBSites()[j]);
				rbcoords[(k*12)+j][0] = x[i]   + rb[0]; 
				rbcoords[(k*12)+j][1] = x[i+1] + rb[1]; 
				rbcoords[(k*12)+j][2] = x[i+2] + rb[2]; 
			}
			k++;
		}
		
		return rbcoords;
		
	}

	private void RMDRVT(int nmol, double[] x, double[][] r, double[][] DRX, double[][] DRY, double[][] DRZ, String ort, boolean g, Rotation rotation) {
		int offset = 3*nmol;
		double[] ri = new double[3];
		double[] p = new double[3];
		
		double[][] RM = new double[3][3];
		double[][] DRMX = new double[3][3];
		double[][] DRMY = new double[3][3];
		double[][] DRMZ = new double[3][3];
		double[] rb = new double[3];
		
		Mat mat = new Mat();
		
		for(int i = 0; i < nmol; i++) {
			int j = 3*(i+1);
			int k = offset + j;
			
			ri[0] = x[j-3];
			ri[1] = x[j-2];
			ri[2] = x[j-1];
			
			p[0] = x[k-3];
			p[1] = x[k-2];
			p[2] = x[k-1];
			
			rotation.rotDrvt(p, RM, DRMX, DRMY, DRMZ, g);

			for(int l = 0; l < 12; l++) {
				int m = 12*i + l;
				
				rb = mat.vdot(RM, getRBSites()[l]);
				r[m][0] = ri[0] + rb[0];
				r[m][1] = ri[1] + rb[1]; 
				r[m][2] = ri[2] + rb[2]; 
				
				if(g) {
				
					rb = mat.vdot(DRMX,  getRBSites()[l]);
					DRX[m][0] = rb[0];
					DRX[m][1] = rb[1];
					DRX[m][2] = rb[2];
				
					rb = mat.vdot(DRMY,  getRBSites()[l]);
					DRY[m][0] = rb[0];
					DRY[m][1] = rb[1];
					DRY[m][2] = rb[2];
				
					rb = mat.vdot(DRMZ,  getRBSites()[l]);
					DRZ[m][0] = rb[0];
					DRZ[m][1] = rb[1];
					DRZ[m][2] = rb[2];
				}
			}
		}
	}

	public double energy(double[] x, Rotation rotation) {
		
		int nmol = (x.length/2)/3;
		double[][] r   = new double[12*nmol][3];
		double[][] DRX = new double[12*nmol][3];
		double[][] DRY = new double[12*nmol][3];
		double[][] DRZ = new double[12*nmol][3];
		String ort = getOrientation();
		
		RMDRVT(nmol, x, r, DRX, DRY, DRZ, ort, false, rotation);
		
		double energy = 0.0;
		
		double[] Rij = new double[3];
		double rij = 0.0;
		double r_1 = 0.0;
		double r_2 = 0.0;
		double r_6 = 0.0;
		
		int o = 0;
		int p = 0;
		
		for(int I = 0; I < nmol-1; I++) {			
			for(int J = I+1; J < nmol; J++) {
				for(int i = 0; i < 12; i++) {
					o 	= (12*I + i);
					
					for(int j = 0; j < 12; j++) {
						p = (12*J + j);
												
						Rij[0] = r[o][0] - r[p][0];
						Rij[1] = r[o][1] - r[p][1];
						Rij[2] = r[o][2] - r[p][2];
						
						rij = Math.sqrt((Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2]));
						r_1 = 1.0/rij;
						r_2 = r_1*r_1;
						r_6 = r_2*r_2*r_2;
						
						if(i < 6 && j < 6) {
							energy += B[0]*Math.exp(-C[0]*rij) - A[0]*r_6 + coulomb*(q[0]*q[0])*r_1;
							//System.out.println(I + " " + J + " " + i + " " + j + ", C-C, " + energy);
						} else if(i >= 6 && j >= 6) {
							energy += B[2]*Math.exp(-C[2]*rij) - A[2]*r_6 + coulomb*(q[1]*q[1])*r_1;
							//System.out.println(I + " " + J + " " + i + " " + j + ", H-H, " + energy);
						} else {
							energy += B[1]*Math.exp(-C[1]*rij) - A[1]*r_6 + coulomb*(q[0]*q[1])*r_1;
							//System.out.println(I + " " + J + " " + i + " " + j + ", C-H, " + energy);
						}
						
					}
				}
			}
		}
		
		return energy;
	}

	public double[] grad(double[] x, Rotation rotation) {
		
		int nmol = (x.length/2)/3;
		double[][] r   = new double[12*nmol][3];
		double[][] DRX = new double[12*nmol][3];
		double[][] DRY = new double[12*nmol][3];
		double[][] DRZ = new double[12*nmol][3];
		
		String ort = getOrientation();
		
		RMDRVT(nmol, x, r, DRX, DRY, DRZ, ort, true,rotation);
		
		double[] g = new double[x.length];
		
		for(int i = 0; i < g.length; i++) {
			g[i] = 0;
		}
		
		double[] Rij = new double[3];
		
		double rij 	= 0.0;
		double r_1 	= 0.0;
		double r_2 	= 0.0;
		double r_7 	= 0.0;
		
		double dUdr = 0.0;
		
		int offset = 3*nmol;
		int o = 0;
		int p = 0;
		
		for(int I = 0; I < nmol-1; I++) {
			int k = 3*(I+1) - 1; //Index for translational coordinate of atom i.
			int l = (offset + k);//Index for rotational coordinate of atom i.
			
			for(int J = I+1; J < nmol; J++) {
				int m = 3*(J+1) - 1; //Index for translational coordinate of atom j.
				int n = (offset + m);//Index for rotational coordinate of atom j.
				
				for(int i = 0; i < 12; i++) {
					o 	= (12*I + i);
					
					for(int j = 0; j < 12; j++) {
						p = (12*J + j);
						
						Rij[0] = r[o][0] - r[p][0];
						Rij[1] = r[o][1] - r[p][1];
						Rij[2] = r[o][2] - r[p][2];
						
						rij 	= Math.sqrt((Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2]));
						r_1 	= 1.0/rij;
						r_2 	= r_1*r_1;
						r_7 	= r_2*r_2*r_2*r_1;
						
						if(i < 6 && j < 6) {
							dUdr = -r_1*(-B[0]*C[0]*Math.exp(-C[0]*rij) + 6.0*A[0]*r_7 - coulomb*(q[0]*q[0])*r_2);
						} else if(i >= 6 && j >= 6) {
							dUdr = -r_1*(-B[2]*C[2]*Math.exp(-C[2]*rij) + 6.0*A[2]*r_7 - coulomb*(q[1]*q[1])*r_2);
						} else {
							dUdr = -r_1*(-B[1]*C[1]*Math.exp(-C[1]*rij) + 6.0*A[1]*r_7 - coulomb*(q[0]*q[1])*r_2);
						}
						
						g[k-2] 	-= dUdr*Rij[0];
						g[k-1] 	-= dUdr*Rij[1];
						g[k] 	-= dUdr*Rij[2];
						
						g[m-2] 	+= dUdr*Rij[0];
						g[m-1] 	+= dUdr*Rij[1];
						g[m] 	+= dUdr*Rij[2];
						
						
						g[l-2] 	-= dUdr*(Rij[0]*DRX[o][0] + Rij[1]*DRX[o][1] + Rij[2]*DRX[o][2]);
						g[l-1] 	-= dUdr*(Rij[0]*DRY[o][0] + Rij[1]*DRY[o][1] + Rij[2]*DRY[o][2]);
						g[l] 	-= dUdr*(Rij[0]*DRZ[o][0] + Rij[1]*DRZ[o][1] + Rij[2]*DRZ[o][2]);
						
						g[n-2] 	+= dUdr*(Rij[0]*DRX[p][0] + Rij[1]*DRX[p][1] + Rij[2]*DRX[p][2]);
						g[n-1] 	+= dUdr*(Rij[0]*DRY[p][0] + Rij[1]*DRY[p][1] + Rij[2]*DRY[p][2]);
						g[n] 	+= dUdr*(Rij[0]*DRZ[p][0] + Rij[1]*DRZ[p][1] + Rij[2]*DRZ[p][2]);
						
					}
				}
			}
		}
		
		return g;
	}
	
	public static void main(String[] args) {
		 
		Benzene b = new Benzene();
		Rotation sxna = new SXNA();
		 
		double[] x = {-0.0804505953, 1.3318831021, 0.3396285110, 
		  	    	   0.0804505953,-2.3318831021,-0.3396285110,
		  	    	   0.8322660626,-0.7372539898, 0.1893163762,
		  	    	   1.8960194774, 1.3595002754, 0.1897039342};
		 
		 System.out.println(b.energy(x, sxna));
		 b.grad(x, sxna);
		 
		 b.viewBZ(x, sxna);
	}

}
