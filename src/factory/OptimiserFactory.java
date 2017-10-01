package factory;

import opt.BFGS;
import opt.ConjGrad;
import opt.LBFGS;
import opt.Optimiser;
import opt.StDes;

public class OptimiserFactory {
	
	public Optimiser chooseOpt(String optimiser) {
		
		
		if(optimiser == "SD") {
			System.out.println("> With Steepest Descent Algorithm");
			return new StDes();
			
		} else if(optimiser == "CG") {
			System.out.println("> With Conjugate Gradient Algorithm");
			return new ConjGrad();
			
		} else if(optimiser == "BFGS") {
			System.out.println("> With BFGS Algorithm");
			return new BFGS();
			
		} else {	
			System.out.println("> With LBFGS Algorithm");
			return new LBFGS();
		}
	}

}
