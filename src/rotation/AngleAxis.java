package rotation;

public class AngleAxis extends Rotation {
	
	public void rotDrvt(double[] p, double[][] RM, double[][] RX, double[][] RY, double[][] RZ, boolean gtest) {
		
		double theta = Math.sqrt((p[0]*p[0]) + (p[1]*p[1]) + (p[2]*p[2])); //Define rotation angle.
		double theta3 = theta*theta*theta;
		double ct = Math.cos(theta);
		double st = Math.sin(theta);
		
		double[] pn = new double[3];	
		
		double[][] pSkew = new double[3][3];
		double[][] pSkew2 = new double[3][3];
		double[][] dPSkew = new double[3][3];
		
		double[][] pSq = new double[3][3];
		double[][] pKp = new double[3][3];
		double[][] ppK = new double[3][3];
		
		Mat mat = new Mat();
		
//---------------------------------------------------------------------------------------------------------------------------------		
//  	Define rotation axis.
//---------------------------------------------------------------------------------------------------------------------------------
		pn[0] = p[0]/theta; pn[1] = p[1]/theta; pn[2] = p[2]/theta;

		
//=================================================================================================================================	
//---------------------------------------------------------------------------------------------------------------------------------		
//  	Compute rotation matrix (RM)
//---------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------
		if(theta == 0) {
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 3; j++) {
					if(i == j) {
						RM[i][j] = 1.0;
					} else {
						RM[i][j] = 0.0;
					}
				}
			}
		} else {
//---------------------------------------------------------------------------------------------------------------------------------
//	  	Find the skew-symmetric matrix of the rotation vector.
//---------------------------------------------------------------------------------------------------------------------------------
			for(int i = 0; i < 3; i++) {
				pSkew[i][i] = 0.0;
			}
			pSkew[0][0] = 0; pSkew[1][1] = 0; pSkew[2][2] = 0;
			pSkew[0][1] = -pn[2]; pSkew[1][0] = pn[2];
			pSkew[0][2] = pn[1]; pSkew[2][0] = -pn[1];
			pSkew[1][2] = -pn[0]; pSkew[2][1] = pn[0];
			
			pSkew2 = mat.dot(pSkew, pSkew);
			
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 3; j++) {
					if(i == j) {
						RM[i][j] = 1.0 + (1.0-ct)*pSkew2[i][j] + st*pSkew[i][j];
					} else {
						RM[i][j] = (1.0-ct)*pSkew2[i][j] + st*pSkew[i][j];
					}
				}
			}
		}
//---------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------
		
	
		if(gtest) {
//=================================================================================================================================			
//---------------------------------------------------------------------------------------------------------------------------------		
//		Derivative of RM w.r.t. x-component of rotation vector.		
//---------------------------------------------------------------------------------------------------------------------------------		
//---------------------------------------------------------------------------------------------------------------------------------		
	
//	Find the derivative of the skew-symmetric matrix of the rotation vector w.r.t. the x-component of the rotation vector.
			dPSkew[0][0] = 0; dPSkew[1][1] = 0; dPSkew[2][2] = 0;
			dPSkew[0][1] = ((p[0]*p[2])/theta3); dPSkew[1][0] = -dPSkew[0][1]; 
			dPSkew[0][2] = -((p[0]*p[1])/theta3); dPSkew[2][0] = -dPSkew[0][2]; 
			dPSkew[1][2] = -((1.0/theta) - ((p[0]*p[0])/theta3)); dPSkew[2][1] = -dPSkew[1][2]; 
			
			pSq = mat.dot(pSkew, pSkew);
			pKp = mat.dot(dPSkew, pSkew);
			ppK = mat.dot(pSkew, dPSkew);
			
//  Compute the derivative of the rotation matrix using Rodrigue's rotation formula.
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 3; j++) {
					RX[i][j] = pn[0]*st*pSq[i][j] + 
							   (1 - ct)*(pKp[i][j] + ppK[i][j]) + 
							   pn[0]*ct*pSkew[i][j] + 
							   st*dPSkew[i][j];
				}
			}
//---------------------------------------------------------------------------------------------------------------------------------		
//---------------------------------------------------------------------------------------------------------------------------------	
			
			
//=================================================================================================================================		
//---------------------------------------------------------------------------------------------------------------------------------		
//				Derivative of RM w.r.t. y-component of rotation vector.		
//---------------------------------------------------------------------------------------------------------------------------------		
//---------------------------------------------------------------------------------------------------------------------------------	
	
//  Find the derivative of the skew-symmetric matrix of the rotation vector w.r.t. the y-component of the rotation vector.
			dPSkew[0][0] = 0; dPSkew[1][1] = 0; dPSkew[2][2] = 0;	
			dPSkew[0][1] = ((p[2]*p[1])/theta3); dPSkew[1][0] = -dPSkew[0][1];
			dPSkew[0][2] = (1.0/theta) - ((p[1]*p[1])/theta3); dPSkew[2][0] = -dPSkew[0][2];	 
			dPSkew[1][2] = ((p[1]*p[0])/theta3); dPSkew[2][1] = -dPSkew[1][2]; 
			
			
			pSq = mat.dot(pSkew, pSkew);
			pKp = mat.dot(dPSkew, pSkew);
			ppK = mat.dot(pSkew, dPSkew);
			
//  Compute the derivative of the rotation matrix using Rodrigue's rotation formula.
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 3; j++) {
					RY[i][j] = pn[1]*st*pSq[i][j] + 
							   (1 - ct)*(pKp[i][j] + ppK[i][j]) + 
							   pn[1]*ct*pSkew[i][j] + 
							   st*dPSkew[i][j];
				}
			}
//---------------------------------------------------------------------------------------------------------------------------------		
//---------------------------------------------------------------------------------------------------------------------------------	
			
					
//=================================================================================================================================			
//---------------------------------------------------------------------------------------------------------------------------------		
//						Derivative of RM w.r.t. z-component of rotation vector.		
//---------------------------------------------------------------------------------------------------------------------------------		
//---------------------------------------------------------------------------------------------------------------------------------	
			
//  Find the derivative of the skew-symmetric matrix of the rotation vector w.r.t. the y-component of the rotation vector.
			dPSkew[0][0] = 0; dPSkew[1][1] = 0; dPSkew[2][2] = 0;
			dPSkew[0][1] = -((1.0/theta) - (p[2]*p[2])/theta3); dPSkew[1][0] = -dPSkew[0][1];  
			dPSkew[0][2] = -((p[1]*p[2])/theta3); dPSkew[2][0] = -dPSkew[0][2];	
		    dPSkew[1][2] = ((p[2]*p[0])/theta3); dPSkew[2][1] = -dPSkew[1][2]; 
			
			pSq = mat.dot(pSkew, pSkew);
			pKp = mat.dot(dPSkew, pSkew);
			ppK = mat.dot(pSkew, dPSkew);
			
//  Compute the derivative of the rotation matrix using Rodrigue's rotation formula.
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 3; j++) {
					RZ[i][j] = pn[2]*st*pSq[i][j] + 
							   (1 - ct)*(pKp[i][j] + ppK[i][j]) + 
							   pn[2]*ct*pSkew[i][j] + 
							   st*dPSkew[i][j];
				}
			}
		}
//---------------------------------------------------------------------------------------------------------------------------------		
//---------------------------------------------------------------------------------------------------------------------------------	
	}

	
	public static void main (String[] args) {
		
		AngleAxis s = new AngleAxis();
		double[] p = {1,0,0};
		double[][] RM = new double[3][3];
		double[][] DRMX = new double[3][3];
		double[][] DRMY = new double[3][3];
		double[][] DRMZ = new double[3][3];
		
		s.rotDrvt(p, RM, DRMX, DRMY, DRMZ,true);
		
		System.out.println(RM[0][0]);
		
	}
}
