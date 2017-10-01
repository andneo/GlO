package rotation;

public class Mat {
	
	public double[][] dot(double[][] a, double[][] b) {
		
		double[][] c = new double[a.length][b[0].length];
		
		for(int i = 0; i < a.length; i++) {
			for(int j = 0; j < b[0].length; j++) {
				for(int k = 0; k < a[0].length; k++) {
					c[i][j] += a[i][k]*b[k][j];
				}
			}
		}
		return c;
	}
	
	public double[] vdot(double[][] a, double[] b) {
		
		double[] c = new double[a.length];
		
		for(int i = 0; i < a.length; i++) {
			for(int j = 0; j < b.length; j++) {
				c[i] += a[i][j]*b[j];
			}
		}
		return c;
	}
	
	public double[][] transpose(double[] a) {
		double[][] b = new double[a.length][1];
		
		for(int i = 0; i < a.length; i++) {
			b[i][0] = a[i];
		}
		
		return b;
	}
}
