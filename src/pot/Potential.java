package pot;

import rotation.Rotation;

/**
 * 
 * @author neo
 * Abstract method for systems of different potential energy models.
 * All systems have an energy based on their geometry, this is given by the energy method.
 * The force acting on all atoms/molecules can be represented by the gradient vector -dU/dX.
 	> U corresponds to the potential energy.
 	> X corresponds to a vector that defines the geometry of the system (e.g. translational or rotational coordinates) and so is a vector.
   	> As a result the grad method returns an array that corresponds to this vector.
 */

public abstract class Potential {

	public abstract double energy(double[] x, Rotation rotation);
	public abstract double[] grad(double[] x, Rotation rotation);

}
