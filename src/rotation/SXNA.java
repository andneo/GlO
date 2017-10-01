package rotation;

import java.util.Random;

public class SXNA extends Rotation {

	public void rotDrvt(double[] p, double[][] RM, double[][] RX, double[][] RY, double[][] RZ, boolean gtest) {
		/*
		 */
		double[] pn = new double[3];
		double[] q = new double[4];
		double theta = Math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
		
		double[][] DRDQ0 = new double[3][3];
		double[][] DRDQ1 = new double[3][3];
		double[][] DRDQ2 = new double[3][3];
		double[][] DRDQ3 = new double[3][3];
		
		double[] DQDX = new double[4];
		double[] DQDY = new double[4];
		double[] DQDZ = new double[4];
		
		pn[0] = p[0]/theta;
		pn[1] = p[1]/theta;
		pn[2] = p[2]/theta;
		
		q[0] = Math.cos(theta/2.0);
		q[1] = (Math.sin(theta/2.0))*pn[0];
		q[2] = (Math.sin(theta/2.0))*pn[1];
		q[3] = (Math.sin(theta/2.0))*pn[2];
		
		RM[0][0] = (q[0]*q[0]) + (q[1]*q[1]) - (q[2]*q[2]) - (q[3]*q[3]);
		RM[0][1] = 2.0*(q[1]*q[2] - q[0]*q[3]);
		RM[0][2] = 2.0*(q[1]*q[3] + q[0]*q[2]);

		RM[1][0] = 2.0*(q[1]*q[2] + q[0]*q[3]);
		RM[1][1] = (q[0]*q[0]) - (q[1]*q[1]) + (q[2]*q[2]) - (q[3]*q[3]);
		RM[1][2] = 2.0*(q[2]*q[3] - q[0]*q[1]);

	    RM[2][0] = 2.0*(q[1]*q[3] - q[0]*q[2]);
		RM[2][1] = 2.0*(q[2]*q[3] + q[0]*q[1]);
		RM[2][2] = (q[0]*q[0]) - (q[1]*q[1]) - (q[2]*q[2]) + (q[3]*q[3]);
		
		if(gtest) {

	        DRDQ0[0][0] = 2.0*q[0];
	        DRDQ0[0][1] = -2.0*q[3];
	        DRDQ0[0][2] = 2.0*q[2];
	
	        DRDQ0[1][0] = 2.0*q[3];
	        DRDQ0[1][1] = 2.0*q[0];
	        DRDQ0[1][2] = -2.0*q[1];
	
	        DRDQ0[2][0] = -2.0*q[2];
	        DRDQ0[2][1] = 2.0*q[1];
	        DRDQ0[2][2] = 2.0*q[0];
	        
	        
	        DRDQ1[0][0] = 2.0*q[1];
	        DRDQ1[0][1] = 2.0*q[2];
	        DRDQ1[0][2] = 2.0*q[3];
	
	        DRDQ1[1][0] = 2.0*q[2];
	        DRDQ1[1][1] = -2.0*q[1];
	        DRDQ1[1][2] = -2.0*q[0];
	
	        DRDQ1[2][0] = 2.0*q[3];
	        DRDQ1[2][1] = 2.0*q[0];
	        DRDQ1[2][2] = -2.0*q[1];
	        
	        
	        DRDQ2[0][0] = -2.0*q[2];
	        DRDQ2[0][1] = 2.0*q[1];
	        DRDQ2[0][2] = 2.0*q[0];
	
	        DRDQ2[1][0] = 2.0*q[1];
	        DRDQ2[1][1] = 2.0*q[2];
	        DRDQ2[1][2] = 2.0*q[3];
	
	        DRDQ2[2][0] = -2.0*q[0];
	        DRDQ2[2][1] = 2.0*q[3];
	        DRDQ2[2][2] = -2.0*q[2];
	        
	        
	        DRDQ3[0][0] = -2.0*q[3];
	        DRDQ3[0][1] = -2.0*q[0];
	        DRDQ3[0][2] = 2.0*q[1];
	
	        DRDQ3[1][0] = 2.0*q[0];
	        DRDQ3[1][1] = -2.0*q[3];
	        DRDQ3[1][2] = 2.0*q[2];
	
	        DRDQ3[2][0] = 2.0*q[1];
	        DRDQ3[2][1] = 2.0*q[2];
	        DRDQ3[2][2] = 2.0*q[3];
	        
	        double thlf = theta/2.0;
	        double tsq = theta*theta;
	        double tcu = theta*theta*theta;
	
	        double a = (Math.sin(thlf))/(2.0*theta);
	        double b = ((Math.cos(thlf)/(2.0*tsq)) - (Math.sin(thlf)/tcu)); 
	        
	        DQDX[0] = -(p[0]*a);
	        DQDX[1] = 2.0*a + (p[0]*p[0])*b;
	        DQDX[2] = (p[0]*p[1])*b;
	        DQDX[3] = (p[0]*p[2])*b;
	        
	        DQDY[0] = -(p[1]*a);
	        DQDY[1] = (p[1]*p[0])*b;
	        DQDY[2] = 2.0*a + (p[1]*p[1])*b;
	        DQDY[3] = (p[1]*p[2])*b;
	        
	        DQDZ[0] = -(p[2]*a);
	        DQDZ[1] = (p[2]*p[0])*b;
	        DQDZ[2] = (p[2]*p[1])*b;
	        DQDZ[3] = 2.0*a + (p[2]*p[2])*b;
	       
	        
	        for(int i = 0; i < 3; i++) {
	        	for(int j = 0; j < 3; j++) {
	        		RX[i][j] = DRDQ0[i][j]*DQDX[0] + DRDQ1[i][j]*DQDX[1] + DRDQ2[i][j]*DQDX[2] + DRDQ3[i][j]*DQDX[3];
	        		RY[i][j] = DRDQ0[i][j]*DQDY[0] + DRDQ1[i][j]*DQDY[1] + DRDQ2[i][j]*DQDY[2] + DRDQ3[i][j]*DQDY[3];
	        		RZ[i][j] = DRDQ0[i][j]*DQDZ[0] + DRDQ1[i][j]*DQDZ[1] + DRDQ2[i][j]*DQDZ[2] + DRDQ3[i][j]*DQDZ[3];
	        	}
	        }
		}
	}
	
	public static void main(String[] args) {
		SXNA s = new SXNA();
		double[] p = {1,0,0};
		double[][] RM = new double[3][3];
		double[][] DRMX = new double[3][3];
		double[][] DRMY = new double[3][3];
		double[][] DRMZ = new double[3][3];
		
		s.rotDrvt(p, RM, DRMX, DRMY, DRMZ,true);
		
		//System.out.println(RM[0][0]);
		
		Random r = new Random(3);
		double[] v = new double[3];
		
		
		for(int i = 0; i < 1000; i++) {
			double azi = 2.0*Math.PI*r.nextDouble();
			double zen = Math.acos(2.0*r.nextDouble() - 1.0);
			
			//System.out.println(Math.asin(Math.sqrt(u)));
			//System.out.println(Math.acos(2.0*u - 1.0) + "\n");
			
			v[0] = Math.sin(zen)*Math.cos(azi);
			v[1] = Math.sin(zen)*Math.sin(azi);
			v[2] = Math.cos(zen);
			
			double norm = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
			
			System.out.println(v[0] + " " + v[1] + " " + v[2]);
			System.out.println(norm+ "\n");
		
		}
	}
}
