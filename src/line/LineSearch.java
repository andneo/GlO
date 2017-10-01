package line;

import pot.Potential;
import rotation.Rotation;

/**
 * 
 * @author neo
 * Simply abstract method for inexact line search methods.
 * See each method for further information.
 */

public abstract class LineSearch {
	
	public abstract double lS(double c1, double c2, double[] g, double[] x, double[] d, Potential pot, Rotation rotation);

}
