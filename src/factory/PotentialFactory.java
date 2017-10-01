package factory;

import pot.Benzene;
import pot.LJ;
import pot.Potential;
import pot.TIP4P;

public class PotentialFactory {
	
	public Potential choosePot(String potential) {
		
		if(potential == "LJ") {
			System.out.println("> Optimising Lennard-Jones System");
			return new LJ();
		} else if(potential == "Benzene") {
			System.out.println("> Optimising Rigid Benzene System");
			return new Benzene();	
		} else {
			System.out.println("> Optimising TIP4P System");
			return new TIP4P();
		}
	}

}
