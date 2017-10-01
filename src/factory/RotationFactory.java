package factory;

import rotation.AngleAxis;
import rotation.Rotation;
import rotation.SXNA;

public class RotationFactory {
	
public Rotation chooseRotation(String rotation) {
		
		
		if(rotation == "SXNA") {
			System.out.println("> Orientation defined by SXNA");
			return new SXNA();
			
		} else {	
			System.out.println("> Orientation defined by Angle Axis");
			return new AngleAxis();
		}
	}

}
