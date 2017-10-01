package gop;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import factory.LineSearchFactory;
import factory.OptimiserFactory;
import factory.PotentialFactory;
import factory.RotationFactory;
import line.LineSearch;
import opt.Optimiser;
import pot.Potential;
import rotation.Rotation;


public class Hop {//java.util.Observable {
	
	public static int iter;
	public static int ite;
	public static boolean fail;
	
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*
 SOFTWARE TO PERFORM GLOBAL OPTIMISATION 
 
 - Currently only basin-hopping is implemented.
 
 - Two chemical systems have been implemented:
 	1. Lennard-Jones clusters with a sigma of 1
 	2. TIP4P model of water
 	
 - Three local optimisation methods are currently implemented:
 	1. Steepest Descent
 	2. Conjugate Gradient
 	3. L-BFGS
 	
 - Each of these methods can use one of three line search methods:
 	1. Backtracking
 	2. Weak Wolfe
 	3. Wolfe
 	
 - For TIP4P cluster the orientaion can be represented by:
 	1. Angle axis coordinates using Rodrigues' rotation formula.
 	2. Angle axis coordinates using exponential map of quaternions.
 */
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private int n;
	private int nAtoms;
	private double[] x;
	private double[] xTemp;
	private double e;
	private double eTemp;
	private double eCurr;
	private int seed;
	private double c;
	private double tStep;
	private double rStep;
	private int stepLow;
	private String LOpt;
	private String potential;
	private String lineSearch;
	private String orientation;
	private boolean random;
	private Random r;
	private double time;
	
	public void setN(int newN) {n = newN;}	
	public int getN() {return n;}
	
	public void setAtoms(int newAtoms) {nAtoms = newAtoms;}
	public int getAtoms() {return nAtoms;}
	
	public void setSeed(int newSeed) {seed = newSeed;}
	public int getSeed() {return seed;}
	
	public void setC(double newC) {c = newC;}
	public double getC() {return c;}
	
	public void setTStep(double newTStep) {tStep = newTStep;}	
	public double getTStep() {return tStep;}
	
	public void setRStep(double newRStep) {
		
		if(newRStep > 2.0*Math.PI) {
			rStep = 2.0*Math.PI;
		} else {
			rStep = newRStep; 
		}
	}
		
	public double getRStep() {return rStep;}
	
	public void setStepLow(int newStepLow) {stepLow = newStepLow;}	
	public int getStepLow() {return stepLow;}
	
	public void setX() {
		if(getPotential() == "LJ") {
			x = new double[3*getAtoms()];
		} else {
			x = new double[6*getAtoms()];
		}
	}
	public void setNewX(double[] newX) {
		for(int j = 0; j < newX.length; j++) {
			x[j] = newX[j];
		}
	}
	public double[] getX() {return x;}
	
	public void setXTemp() {
		if(getPotential() == "LJ") {
			xTemp = new double[3*getAtoms()];
		} else {
			xTemp = new double[6*getAtoms()];
		}
	}
	public void setNewXTemp(double[] newXTemp) {
		for(int j = 0; j < newXTemp.length;j++) {
			xTemp[j] = newXTemp[j];
		}
	}
	public double[] getXTemp() {return xTemp;}
	
	public void setE(double newE) {e = newE;}	
	public double getE() {return e;}
	
	public void setETemp(double newETemp) {eTemp = newETemp;}	
	public double getETemp() {return eTemp;}
	
	public void setECurr(double newECurr) {eCurr = newECurr;}	
	public double getECurr() {return eCurr;}
	
	public void setOptimiser(String newOpt) {LOpt = newOpt;}	
	public String getOptimiser() {return LOpt;}
	
	public void setLS(String newLS) {lineSearch = newLS;}	
	public String getLS() {return lineSearch;}
	
	public void setPotential(String newPot) {potential = newPot;}	
	public String getPotential() {return potential;}
	
	public String getOrient() { return orientation; }
	public void setOrient(String o) { orientation = o;}
	
	public void setRandom(boolean newR) {random = newR;}
	public boolean getR() {return random;}
	
	public void setTime(double d) {time = d;}	
	public double getTime() {return time;}
	
	public void setRand() {r = new Random(getSeed());}
	public Random getRand() {return r;}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------                                                       BASIN-HOPPING METHOD BEGINS HERE                                                                              ------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public double gOptimise(int n, int nAtoms, double c, int seed, String LOpt, String lineSearch, String potential, String orientation, boolean random, double tStep, double rStep) {
		setN(n); 				 //Number of Monte Carlo steps.
		setPotential(potential); //System being modelled (LJ or TIP4P)
		setAtoms(nAtoms); 		 //Size of the system (no. atoms/molecules)
		setC(c); 				 //Size of container
		setSeed(seed); 			 //Seed for pseudo-random number generator.
		setOptimiser(LOpt); 	 //Optimisation method used. 
		setLS(lineSearch); 		 //Line search method used.
		setRandom(random);		 //True for completely random steps
		setTStep(tStep);         //Set translational step size.
		setRStep(rStep);         //Set rotational step size.
		setRand();				 //Initialise pseudo-random number generator.
		setOrient(orientation);
		
		setX(); 				 // Coordinates of current Monte Carlo step.
		setXTemp(); 			 // Coordinates of last successful Monte Carlo step.
		setE(0); 				 // Energy of current lowest minimum found.
		setETemp(0); 			 // Energy of last successful Monte Carlo step.
		setECurr(0); 			 // Energy of current Monte Carlo step.
		
		/*
		 * Create factories and corresponding objects.
		 * This enables the system being modelled and the optimisation method used to be easily interchanged. 
		 */

		System.out.println("> Generating Factories......");
		PotentialFactory potFac = new PotentialFactory();
		Potential pot = potFac.choosePot(getPotential());
		
		OptimiserFactory optFac = new OptimiserFactory();
		Optimiser opt = optFac.chooseOpt(getOptimiser());
		
		LineSearchFactory lSF = new LineSearchFactory();
		LineSearch lS = lSF.chooseLS(getLS());
		
		RotationFactory rf = new RotationFactory();
		Rotation rotation = rf.chooseRotation(getOrient());
		
		Writer wo = new Writer();
		Centre centre = new Centre();
		Rand rand = new Rand();
		
		System.out.println();
		
		//Token array used when moving centre of geometry to the origin
		double[] xCentre = new double[getX().length];

		boolean t1 = false;
		boolean newMin = true;
		int failure = 0; // Counter for the number of steps a new 'global' minimum hasn't been found

		int i1 = 0; // Counter for the total number of Monte Carlo steps.
		
		//List of the lowest minima found.
		ArrayList<Minimum> lowest = new ArrayList<Minimum>(); 
		
		double time = 0;		
		time = System.currentTimeMillis();
		
		System.out.println("> No. Atoms: " + getAtoms());
		System.out.println("> T Step Size: " + getTStep());
		System.out.println("> R Step Size: " + getRStep());
		System.out.println("> Cell Size: " + getC());
		System.out.println("> Seed: " + getSeed() + "\n");
		
		fail = false;
		
		wo.start(getOptimiser(), getAtoms(), getN(), getSeed(), getC(), getTStep(), getRStep(), getLS(), getOrient());
		
		Random r1 = new Random(getSeed());
		
		for(int i = 0; i <= getN(); i++) {
			if(getR()) {
				t1 = false;
			} else if(failure == getAtoms()*500) {
				t1 = false;
				failure = 0;
				wo.failure(getOptimiser(), getAtoms(), getN(), getSeed());
			}
			
			fail = false;
			/*
			 * Choose a set of initial coordinates the atoms.
			 * A set of pseudo-random doubles between -1 and 1 are chosen.
			 * If random basin-hopping is selected then a set of new random coordinates are chosen at each step.
			 */
			while(t1 == false) { 
				setNewX(rand.setRandX(getX(), getC(), getAtoms(), getPotential(), getRand()));
				fail = false;
				
				// If one (or more) atoms occupy the same position a new set of coordinates must be chosen.
				if(Double.isNaN(getETemp())) {
					t1 = false;
				} else {
					t1 = true;
				}		
			}
			
		   /* 
		    * Randomly perturb the system.
			* A set of pseudo-random doubles between -1 and 1 are chosen.
			* This double is then multiplied by the step size.
			*/
			if(i > 0 && getR() == false) {
				setNewX(rand.setRandX2(getX(), getTStep(), getRStep(), getPotential(), getRand()));
			}
			
			//Optimise the new perturbed geometry
			setNewX(opt.optimise(getX(), lS, pot, rotation));
			//Move the centre of geometry to the origin
			xCentre = centre.centre(getX(), getAtoms());
			setNewX(xCentre);
			
			//Calculate the energy of the new geometry.
			setECurr(pot.energy(getX(), rotation));
			 
			double a = r1.nextDouble(); 
			
			/*
			 * 1.If the current energy is lower than eTemp the step is accepted.
			 *   Or if the Metropolis condition is fulfilled accept the step.
			 * 2.If accepted set eTemp to be the current energy.
			 * 3.Update xTemp with the current coordinates.
			 * 	 This is done in order to be able to reset x when a step fails.
			 */
			
			if(getECurr() < getETemp() || a < Math.exp(-(getECurr() - getETemp())/2.0)) { //1
				
				setETemp(getECurr()); //2
				setNewXTemp(getX()); //3
				

				if(i == 0) {
					setE(getETemp());
				} else {
					if(getETemp() < getE()) {
						setE(getETemp());
						setStepLow(i);
						failure = 0;
					} else {
						failure++;
					}
				}
				
				if(lowest.size() < 11) {
					for(int j = 0; j < lowest.size(); j++) {
						if(Math.abs(lowest.get(j).getE()-getETemp()) < 0.01) {
							newMin = false;
							break;
						} 
					}
					
					if(newMin) {
						double[] temp = new double[getX().length];
						for(int l = 0; l < getX().length; l++) {
							temp[l] = getX()[l];
						}
						Minimum m = new Minimum(temp,getETemp(),i);
						lowest.add(m);
						Collections.sort(lowest);
					}
					newMin = true;
				} else if(lowest.get(10).getE() > getETemp()) {	
					for(int j = 0; j < lowest.size(); j++) {
						if(Math.abs(lowest.get(j).getE()-getETemp()) < 0.01) {
							newMin = false;
							break;
						} 
					}
					if(newMin) {
						lowest.remove(10);
						double[] temp = new double[getX().length];
						for(int l = 0; l < getX().length; l++) {
							temp[l] = getX()[l];
						}
						Minimum m = new Minimum(temp,getETemp(),i);
						lowest.add(m);
						Collections.sort(lowest);
					}
					newMin = true;
				}
				
				//wo.sFilm(getAtoms(), getX(), getPotential(), rotation, i);
				
				wo.acceptedStep(getOptimiser(), getAtoms(), getN(), getSeed(), getC(), 
						        getE(), getETemp(), i1, ite, i, (System.currentTimeMillis() - time));
				
				i1++;
			
			} else {
				//If a step is not accepted increment the fail counter.
				failure++;
				//wo.fFilm(getAtoms(), getX(), getPotential(), rotation, i);
				//Reset the x to the previously accepted step.
				setNewX(getXTemp());
				
				wo.rejectedStep(getOptimiser(), getAtoms(), getN(), getSeed(), getC(), 
						        getE(), getECurr(), i1, ite, i, (System.currentTimeMillis() - time));
			}
			
			//If a step fails write it out in the output.
			if(fail == true) {
				wo.failedStep(getOptimiser(), getAtoms(), getN(), getSeed(),
				              ite, (System.currentTimeMillis() - time));
			}
			
			System.gc();
		}
		
		wo.finish(getOptimiser(), getAtoms(), getN(), getSeed(), iter, 
				  getE(), getStepLow(), (System.currentTimeMillis() - time));
		
		// Write coordinates of lowest 10 structures found to file.
		wo.coords(getOptimiser(), getAtoms(), getN(), getSeed(), lowest, getPotential(), rotation);
		
		setTime((System.currentTimeMillis() - time));
		
		System.out.println("> Lowest Energy Found: " + getE());
    	System.out.println("> Average Number of LOpt Steps: " + (iter/getN()));
    	System.out.println("> Time Taken to Complete: " + (getTime()*0.001) + "s");
		iter = 0;
		return getE();
	
	}
	
	public static void iter(int it) {
		
		iter += it; //Total number of local optimisation steps.
		ite = it;   //Current number of local optimisation steps.
		
	}
	
	public static void main(String[] args) {
		
		//String opt = "SD";
		//String opt = "CG";
		//String opt = "BFGS";
		String opt = "LBFGS";
		
		//String lS = "BackTrack";
		//String lS = "WeakWolfe";
		String lS = "Wolfe";
		
		String pot = "LJ";
		//pot = "TIP";
		
		//String o = "AA";
		String o = "SXNA";
		
		Hop s = new Hop();
		
		//setParameters(int n, int nAtoms, double c, int seed, String LOpt, String lineSearch, String potential, String orientation, boolean random, double tStep, double rStep)
		s.gOptimise(10, 100, 5.0, 7, opt, lS, pot, o, false, 0.5, 1.0);
	}
}