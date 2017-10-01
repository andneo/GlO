package opt;
import line.LineSearch;
import pot.Potential;
import rotation.Rotation;

/**
 * 
 * @author neo
 * Abstract class to define all local optimisation methods implemented.
 * The optimisers take in the vector x as a parameter, this is what is to be optimised.
 	> The vector x corresponds to the geometry of the system that is to be optimised.
 * Different systems have different potential energy functions.
 	> Hence the pot argument ensures the optimiser knows what kind of system is being optimised.
 * The local optimisation methods are gradient based and so the direction of optimisation is defined by the gradient vector.
 	> However the optimal step to take along this direction is unknown.
 	> Various inexact line search methods have been implemented to determine the step length.
 	> The lineSearch argument defines which method is to be used.
 */
public abstract class Optimiser {
	
	public abstract double[] optimise(double[] x, LineSearch lineSearch, Potential pot, Rotation rotation);

}
