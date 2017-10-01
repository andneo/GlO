package factory;

import line.BackTrack;
import line.WeakWolfe;
import line.Wolfe;
import line.LineSearch;

public class LineSearchFactory {
	
	public LineSearch chooseLS(String lS) {
		
		if(lS == "BackTrack") {
			System.out.println("> Using Back Tracking Line Search");
			return new BackTrack();
		} else if(lS == "WeakWolfe") {
			System.out.println("> Using Weak Wolfe Line Search");
			return new WeakWolfe();
		} else {
			System.out.println("> Using Wolfe Line Search");
			return new Wolfe();
		}
	}

}
