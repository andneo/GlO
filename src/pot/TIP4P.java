package pot;

import rotation.Mat;
import rotation.Rotation;

public class TIP4P extends Potential {
	
	private double theta = (Math.PI*104.52)/180.0;
	
	private double[][] rbsites = {{0.0,0.0,0.0},
								  {0.0,Math.sin(0.5*theta)*0.9572,Math.cos(0.5*theta)*0.9572},
								  {0.0,-Math.sin(0.5*theta)*0.9572,Math.cos(0.5*theta)*0.9572},
								  {0.0,0.0,0.15}};
	
	private double[] charge = {0.0,0.52,0.52,-1.04};
	
	private double A = 2510.4E03;
	private double B = 2552.24;
	private double coulomb = 1389.354848;
	
	private double[] mass = {16.0, 1.0, 1.0, 0.0};
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
		
		for(int i = 0; i < 4; i++) {
			cm[0] = cm[0] + mass[i]*rbsites[i][0]; 
			cm[1] = cm[1] + mass[i]*rbsites[i][1]; 
			cm[2] = cm[2] + mass[i]*rbsites[i][2];
			m += mass[i];
		}
		
		for(int i = 0; i < 3; i++) {
			cm[i] = cm[i]/m; 
		}
		
		for(int i = 0; i < 4; i++) {
			rbsites[i][0] = rbsites[i][0] - cm[0];
			rbsites[i][1] = rbsites[i][1] - cm[1];
			rbsites[i][2] = rbsites[i][2] - cm[2];
		}
		
		return rbsites;	
	}
	
	public double[][] viewTIP(double[] x, Rotation rotation) {
		
		double[] p      	= new double[3];
		double[][] RM   	= new double[3][3];
		double[][] DRMX 	= new double[3][3];
		double[][] DRMY 	= new double[3][3];
		double[][] DRMZ 	= new double[3][3];
		double[] rb     	= new double[3];
		double[][] rbcoords = new double[((x.length/2)/3)*3][3];
		int offset 			= x.length/2;
		
		for(int i = 0; i < offset; i+=3) {
			p[0] = x[i+offset];
			p[1] = x[i+1+offset];
			p[2] = x[i+2+offset];
			
			Mat mat = new Mat();
			rotation.rotDrvt(p, RM, DRMX, DRMY, DRMZ, false);
			
			for(int j = 0; j < 3; j++) {
				rb = mat.vdot(RM,  getRBSites()[j]);
				rbcoords[i+j][0] = x[i] + rb[0]; 
				rbcoords[i+j][1] = x[i+1] + rb[1]; 
				rbcoords[i+j][2] = x[i+2] + rb[2]; 
			}
		}
		
		return rbcoords;
		
	}
	
	private void RMDRVT(int nmol, double[] x, double[][] r, double[][] DRX, double[][] DRY, double[][] DRZ, String ort, boolean g, Rotation rotation) {
		int offset 	= 3*nmol;
		double[] ri = new double[3];
		double[] p 	= new double[3];
		
		double[][] RM 	= new double[3][3];
		double[][] DRMX = new double[3][3];
		double[][] DRMY = new double[3][3];
		double[][] DRMZ = new double[3][3];
		double[] rb 	= new double[3];
		
		Mat mat = new Mat();
		
		for(int i = 0; i < nmol; i++) {
			int j = 3*(i+1);
			int k = offset + j;
			
			//Coordinates for center of mass of molecule
			ri[0] = x[j-3]; 
			ri[1] = x[j-2];
			ri[2] = x[j-1];

			//Rotational coordinates of molecule
			p[0] = x[k-3];
			p[1] = x[k-2];
			p[2] = x[k-1];
			
			//Calculate rotation matrix (RM) for molecule
			rotation.rotDrvt(p, RM, DRMX, DRMY, DRMZ, g);

			for(int l = 0; l < 4; l++) {
				int m = 4*i + l;
				
				//Calculate absolute position of sites in molecule
				// r_i = r_I + RM_I*r_i0
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
		double[][] r   = new double[4*nmol][3];
		double[][] DRX = new double[4*nmol][3];
		double[][] DRY = new double[4*nmol][3];
		double[][] DRZ = new double[4*nmol][3];
		String ort = getOrientation();
		
		RMDRVT(nmol, x, r, DRX, DRY, DRZ, ort, false, rotation);
		
		double energy = 0.0;
		
		double[] rss = new double[3];
		double r2 = 0.0;
		double r6 = 0.0;
		double r12 = 0.0;
		
		for(int i = 0; i < nmol; i++) {			
			for(int l = i+1; l < nmol; l++) {
				
				int o = (4*i + 1) - 1;
				int p = (4*l + 1) - 1;
				
				rss[0] = r[o][0] - r[p][0];
				rss[1] = r[o][1] - r[p][1];
				rss[2] = r[o][2] - r[p][2];
				
				r2 = 1.0/(rss[0]*rss[0] + rss[1]*rss[1] + rss[2]*rss[2]);
				r6 = r2*r2*r2;
				r12 = r6*r6;
				
				// Lennard-Jones contribution to the energy
				energy += A*r12 - B*r6;
				
				for(int I = 1; I < 4; I++) {
					o = (4*i + I);
					for(int J = 1; J < 4; J++) {
						p = (4*l + J);
						
						rss[0] = r[o][0] - r[p][0];
						rss[1] = r[o][1] - r[p][1];
						rss[2] = r[o][2] - r[p][2];
						r2 = 1.0/Math.sqrt((rss[0]*rss[0] + rss[1]*rss[1] + rss[2]*rss[2]));
						
						// Coulombic contribution to the energy
						energy += coulomb*charge[I]*charge[J]*r2;
					}
				}			
			}
		}
		
		return energy;
	}

	public double[] grad(double[] x, Rotation rotation) {
		
		int nmol = (x.length/2)/3;
		double[][] r = new double[x.length][3];
		double[][] DRX = new double[4*nmol][3];
		double[][] DRY = new double[4*nmol][3];
		double[][] DRZ = new double[4*nmol][3];
		
		double[][] DVDR = new double[4*nmol][4*nmol];
		
		double[][] DOTIX = new double[4*nmol][4*nmol];
		double[][] DOTIY = new double[4*nmol][4*nmol];
		double[][] DOTIZ = new double[4*nmol][4*nmol];
		
		double[][] DOTJX = new double[4*nmol][4*nmol];
		double[][] DOTJY = new double[4*nmol][4*nmol];
		double[][] DOTJZ = new double[4*nmol][4*nmol];
		
		String ort = getOrientation();
		
		RMDRVT(nmol, x, r, DRX, DRY, DRZ, ort, true,rotation);
		
		double[] rss = new double[3];
		double r2 = 0.0;
		double r6 = 0.0;
		double r12 = 0.0;
		double absr = 0.0;
		
		double[] g = new double[x.length];
		
		for(int i = 0; i < g.length; i++) {
			g[i] = 0;
		}
		
		int offset = 3*nmol;
		
		for(int i = 0; i < nmol; i++) {	
			int j = 3*(i+1) - 1; //Index for translational coordinate of atom i.
			int k = (offset + j);//Index for rotational coordinate of atom i.
			
			for(int l = i+1; l < nmol; l++) {
				int m = 3*(l+1) - 1; //Index for translational coordinate of atom j.
				int n = (offset + m);//Index for rotational coordinate of atom j.
				
				int o = (4*i + 1) - 1; 
				int p = (4*l + 1) - 1; 
				
				rss[0] = r[o][0] - r[p][0];
				rss[1] = r[o][1] - r[p][1];
				rss[2] = r[o][2] - r[p][2];
				
				r2 = 1.0/(rss[0]*rss[0] + rss[1]*rss[1] + rss[2]*rss[2]);
				r6 = r2*r2*r2;
				r12 = r6*r6;
				
				DVDR[o][p] = -6.0*(2.0*A*r12 - B*r6)*r2;
				DVDR[p][o] = DVDR[o][p];
				
				DOTIX[o][p] = rss[0]*DRX[o][0] + rss[1]*DRX[o][1] + rss[2]*DRX[o][2];
				DOTIY[o][p] = rss[0]*DRY[o][0] + rss[1]*DRY[o][1] + rss[2]*DRY[o][2];
				DOTIZ[o][p] = rss[0]*DRZ[o][0] + rss[1]*DRZ[o][1] + rss[2]*DRZ[o][2];
				
				DOTJX[o][p] = rss[0]*DRX[p][0] + rss[1]*DRX[p][1] + rss[2]*DRX[p][2];
				DOTJY[o][p] = rss[0]*DRY[p][0] + rss[1]*DRY[p][1] + rss[2]*DRY[p][2];
				DOTJZ[o][p] = rss[0]*DRZ[p][0] + rss[1]*DRZ[p][1] + rss[2]*DRZ[p][2];
				
				g[j-2]  += DVDR[o][p]*rss[0];
				g[j-1]  += DVDR[o][p]*rss[1];
				g[j]    += DVDR[o][p]*rss[2];
				
				g[m-2]  -= DVDR[o][p]*rss[0];
				g[m-1]  -= DVDR[o][p]*rss[1];
				g[m]    -= DVDR[o][p]*rss[2];
				
				g[k-2]  += DVDR[o][p]*DOTIX[o][p];
				g[k-1]  += DVDR[o][p]*DOTIY[o][p];
				g[k]    += DVDR[o][p]*DOTIZ[o][p];
				
				g[n-2]  -= DVDR[o][p]*DOTJX[o][p];
				g[n-1]  -= DVDR[o][p]*DOTJY[o][p];
				g[n]    -= DVDR[o][p]*DOTJZ[o][p];
				
				
				for(int I = 1; I < 4; I++) {
					o = (4*i + I);
					for(int J = 1; J < 4; J++) {
						p = (4*l + J);
						
						rss[0] = r[o][0] - r[p][0];
						rss[1] = r[o][1] - r[p][1];
						rss[2] = r[o][2] - r[p][2];
						absr = Math.sqrt(rss[0]*rss[0] + rss[1]*rss[1] + rss[2]*rss[2]);
						r2 = 1.0/(rss[0]*rss[0] + rss[1]*rss[1] + rss[2]*rss[2]);
						
						DVDR[o][p] = -coulomb*charge[I]*charge[J]*r2/absr;
						DVDR[p][o] = DVDR[o][p];
						
						DOTIX[o][p] = rss[0]*DRX[o][0] + rss[1]*DRX[o][1] + rss[2]*DRX[o][2];
						DOTIY[o][p] = rss[0]*DRY[o][0] + rss[1]*DRY[o][1] + rss[2]*DRY[o][2];
						DOTIZ[o][p] = rss[0]*DRZ[o][0] + rss[1]*DRZ[o][1] + rss[2]*DRZ[o][2];
						
						DOTJX[o][p] = rss[0]*DRX[p][0] + rss[1]*DRX[p][1] + rss[2]*DRX[p][2];
						DOTJY[o][p] = rss[0]*DRY[p][0] + rss[1]*DRY[p][1] + rss[2]*DRY[p][2];
						DOTJZ[o][p] = rss[0]*DRZ[p][0] + rss[1]*DRZ[p][1] + rss[2]*DRZ[p][2];
						
						g[j-2]  += DVDR[o][p]*rss[0];
						g[j-1]  += DVDR[o][p]*rss[1];
						g[j]    += DVDR[o][p]*rss[2];
						
						g[m-2]  -= DVDR[o][p]*rss[0];
						g[m-1]  -= DVDR[o][p]*rss[1];
						g[m]    -= DVDR[o][p]*rss[2];
						
						g[k-2]  += DVDR[o][p]*DOTIX[o][p];
						g[k-1]  += DVDR[o][p]*DOTIY[o][p];
						g[k]    += DVDR[o][p]*DOTIZ[o][p];
						
						g[n-2]  -= DVDR[o][p]*DOTJX[o][p];
						g[n-1]  -= DVDR[o][p]*DOTJY[o][p];
						g[n]    -= DVDR[o][p]*DOTJZ[o][p];
					}
				}		
			}
		}
		
		return g;
	}
}
