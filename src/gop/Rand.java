package gop;

import java.util.Random;

public class Rand {
	
	/*
	 * Setters for random perturbation of system geometry.
	 * RandX generates a completely random set of new coordinates.
	 * RandX2 randomly perturbs the existing geometry of the system.
	 */
	
	public double[] setRandX(double[] x, double c, int atoms, String potential, Random r) {
		
		if(potential == "LJ") {
			for(int a = 0; a < x.length; a++) {
				x[a] = ((r.nextDouble()*2.0) - 1.0*c);
			}
		} else {
			
			int offset = atoms*3;
			double theta = 0.0;
			double qdot = 0.0;
			double[] q = new double[4];
			
			for(int a = 0; a < offset; a += 3) {
				
				// Generate random translational coordinates
				x[a]   = ((r.nextDouble()*2.0) - 1.0)*c;
				x[a+1] = ((r.nextDouble()*2.0) - 1.0)*c;
				x[a+2] = ((r.nextDouble()*2.0) - 1.0)*c;
				
				// Generate random rotational coordinates from random unit quaternion
				q = randQuat(r);
				
				qdot = Math.sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
				theta = 2.0*Math.acos(q[0]);
				
				// Convert unit quaternion into angle-axis vector 
				if(theta <= 1.0D-12) {
					x[a+offset]   = 0.0;
					x[a+1+offset] = 0.0;
					x[a+2+offset] = 0.0;
				} else {
					x[a+offset]   = theta*q[1]/qdot;
					x[a+1+offset] = theta*q[2]/qdot;
					x[a+2+offset] = theta*q[3]/qdot;
				}
			}
		}
		
		return x;
	}
	
	public double[] setRandX2(double[] x, double tStep, double rStep, String potential, Random r) {
		
		double[] newX = new double[x.length];
		
		if(potential == "LJ") {
			for(int a = 0; a < x.length; a++) {
				newX[a] = x[a] + ((r.nextDouble()*2.0) - 1.0)*tStep;
			}
		} else {
			
			double[] p = new double[3];
			double[] q = new double[4];
			double[] q1 = new double[4];
			double[] q2 = new double[4];
			double theta = 0.0; 
			double qdot = 0.0;
			
			double[] axis = new double[3];
			double angle = 0.0;
			double norm = 0.0;
			double azi = 0.0;
			double zen = 0.0;
			double s = 0.0;
			double z = 0.0;
			
			for(int a = 0; a < x.length/2; a++) {
				newX[a] = x[a] + ((r.nextDouble()*2.0) - 1.0)*tStep;
			}
			for(int i = x.length/2; i < x.length; i+=3) {
				p[0] = x[i];
				p[1] = x[i+1];
				p[2] = x[i+2];
			
				theta = Math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
			
				// Convert current rotational coordinates to unit quaternion
				q1[0] = Math.cos(theta/2.0);
				q1[1] = (Math.sin(theta/2.0))*(p[0]/theta);	
				q1[2] = (Math.sin(theta/2.0))*(p[1]/theta);
				q1[3] = (Math.sin(theta/2.0))*(p[2]/theta);
				
		
				// Choose a random axis of rotation by forming random unit vector.
				azi = 2.0*Math.PI*r.nextDouble();
				zen = Math.acos(2.0*r.nextDouble() - 1.0);
				
				axis[0] = Math.sin(zen)*Math.cos(azi);
				axis[1] = Math.sin(zen)*Math.sin(azi);
				axis[2] = Math.cos(zen);		

				norm = Math.sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);	
				if(norm != 1.0) {
					for(int j = 0; j < 3; j++) {
						axis[j] = axis[j]/norm;
					}
				}
				
				// Choose random angle of rotation in range 0:step
				// Angle chosen from a distribution of (sin(0.5*step))**2
				s = 1.0/(Math.sin(0.5*rStep)*Math.sin(0.5*rStep));
				angle = r.nextDouble()*rStep;
				z = r.nextDouble();
				
				while(z > s*(Math.sin(0.5*angle)*Math.sin(0.5*angle))) {
					angle = r.nextDouble()*rStep;
				}
				
				angle = angle/2.0;
				
				//angle = (getRand().nextDouble()*Math.PI)/2.0;
				
				// Convert random angle-axis vector to unit quaternion
				q2[0] = Math.cos(angle);
				q2[1] = (Math.sin(angle))*(axis[0]);	
				q2[2] = (Math.sin(angle))*(axis[1]);
				q2[3] = (Math.sin(angle))*(axis[2]);
			
				// Combine existing rotational coordinates with new random ones
				// This is evaluating the dot product their corresponding quaternions
				q[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
				q[1] = q1[1]*q2[0] + q1[0]*q2[1] + q1[3]*q2[2] - q1[2]*q2[3];
				q[2] = q1[2]*q2[0] + q1[0]*q2[2] + q1[1]*q2[3] - q1[3]*q2[1];
				q[3] = q1[3]*q2[0] + q1[0]*q2[3] + q1[2]*q2[1] - q1[1]*q2[2];
			
				//Convert combined quaternion into angle-axis vector
				theta = 2.0*Math.acos(q[0]);
				qdot = Math.sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
			
				if(theta <= 1.0D-12) {
					newX[i] = 0.0;
					newX[i+1] = 0.0;
					newX[i+2] = 0.0;
				} else {
					newX[i]   = theta*q[1]/qdot;
					newX[i+1] = theta*q[2]/qdot;
					newX[i+2] = theta*q[3]/qdot;
				}
			}
		}
		
		return newX;
	}
	
	public double[] randQuat(Random r) {
		double[] v = new double[3];
		double qdot = 0.0;
		double[] q = new double[4];	
		double p1  = 0.0;
		double p2  = 0.0;
		double a1  = 0.0;
		double a2  = 0.0;	
		double pi2 = Math.PI*2.0;
		
		// Generate random rotational coordinates from random unit quaternion
		v[0] = r.nextDouble();
		v[1] = r.nextDouble();
		v[2] = r.nextDouble();
		
		p1 = Math.sqrt(1.0 - v[0]);
		p2 = Math.sqrt(v[0]);
		
		a1 = pi2*v[1];
		a2 = pi2*v[2];
		
		q[0] = p1*Math.sin(a1);
		q[1] = p1*Math.cos(a1);
		q[2] = p2*Math.sin(a2);
		q[3] = p2*Math.cos(a2);
		
		qdot = Math.sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
		
		// If not a unit quaternion normalise
		if(qdot != 1.0) {
			q[0] = q[0]/qdot;
			q[1] = q[1]/qdot;
			q[2] = q[2]/qdot;
			q[3] = q[3]/qdot;
		} 
		
		return q;
	}

}
