package gop;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import pot.TIP4P;
import rotation.Rotation;

public class Writer {
	
	private BufferedWriter writer;
	
	public void start(String optimiser, int atoms, int n, int seed, double c, double tStep, double rStep, String linesearch, String orientation) {
		try {
			writer = new BufferedWriter(new FileWriter(optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".txt", false));
			 writer.write("NATOMS: " + String.valueOf(atoms));
			 writer.newLine();
			 writer.write("LOCAL MINIMISER: " + optimiser + " with " + linesearch);
			 writer.newLine();
			 writer.write("CELL SIZE: " + String.valueOf(c));
			 writer.newLine();
			 writer.write("TRANSLATIONAL STEP SIZE: " + String.valueOf(tStep));
			 writer.newLine();
			 writer.write("ROTATIONAL STEP SIZE: " + String.valueOf(rStep));
			 writer.newLine();
			 writer.write("ORIENTATION: " + orientation);
			 writer.newLine();
			 writer.write("SEED: " + String.valueOf(seed));
			 writer.newLine();
			 writer.newLine();
		} catch (IOException e1) {
		}
		{
			try {
				if (writer != null)
					writer.close();
			} catch (IOException e1) {
			}
		}
	}
	
	public void failure(String optimiser, int atoms, int n, int seed) {
		try {
			 writer = new BufferedWriter(new FileWriter(optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".txt", true));
			 writer.write("**************************************************************************");
			 writer.newLine();
			 writer.write("New minimum has not been found for " + String.valueOf(atoms*500)+ " steps");
			 writer.newLine();
			 writer.write("Resetting with new random structure");
			 writer.newLine();
			 writer.write("**************************************************************************");
			 writer.newLine();
			 writer.newLine();
		} catch (IOException e1) {
		}
		{
			try {
				if (writer != null)
					writer.close();
			} catch (IOException e1) {
			}
		}
	}
	
	public void acceptedStep(String optimiser, int atoms, int n, int seed, double c, double e, double eTemp, int i, int ite, int i1, double time) {
		try {
			writer = new BufferedWriter(new FileWriter(optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".txt", true));
			writer.write("ACCEPTED STEP"+ " (" + String.valueOf(i) + ")");
			writer.newLine();
			writer.write("LOpt Steps: " + String.valueOf(ite));
			writer.newLine();
			writer.write("LOWEST ENERGY: " + String.valueOf(e));
			writer.newLine();
			writer.write("CURRENT ENERGY: " + String.valueOf(eTemp)+ " (" + String.valueOf(i1) + ")");
			writer.newLine();
			writer.write("TIME: " + String.valueOf((float)(time*0.001)) + "s");
			writer.newLine();
			writer.newLine();
		} catch (IOException e1) {
		}
		{
			try {
				if (writer != null)
					writer.close();
			} catch (IOException e1) {
			}
		}
	}
	
	public void rejectedStep(String optimiser, int atoms, int n, int seed, double c, double e, double eCurr, int i, int ite, int i1, double time) {
		try {
			writer = new BufferedWriter(new FileWriter(optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".txt", true));
			writer.write("REJECTED STEP"+ " (" + String.valueOf(i) + ")");
			writer.newLine();
			writer.write("LOpt Steps: " + String.valueOf(ite));
			writer.newLine();
			writer.write("LOWEST ENERGY: " + String.valueOf(e));
			writer.newLine();
			writer.write("CURRENT ENERGY: " + String.valueOf(eCurr)+ " (" + String.valueOf(i1) + ")");
			writer.newLine();
			writer.write("TIME: " + String.valueOf((float)(time*0.001)) + "s");
			writer.newLine();
			writer.newLine();
		} catch (IOException e1) {
		}
		{
			try {
				if (writer != null)
					writer.close();
			} catch (IOException e1) {
			}
		}
	}
	
	public void failedStep(String optimiser, int atoms, int n, int seed, int ite, double time) {
		try {
			writer = new BufferedWriter(new FileWriter(optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".txt", true));
			writer.write("FAILED STEP");
			writer.newLine();
			writer.write("LOpt Steps: " + String.valueOf(ite));
			writer.newLine();
			writer.write("TIME: " + String.valueOf((float)(time*0.001)) + "s");
			writer.newLine();
			writer.newLine();
		} catch (IOException e1) {
		}
		{
			try {
				if (writer != null)
					writer.close();
			} catch (IOException e1) {
			}
		}
	}
	
	public void finish(String optimiser, int atoms, int n, int seed, int iter, double e, int stepLow, double time) {
		int minutes = (int)(time*0.001)/60;
		try {
			 writer = new BufferedWriter(new FileWriter(optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".txt", true));
			 writer.write("TIME: " + String.valueOf(minutes) + " minutes and " + String.valueOf((float)(time*0.001) - (minutes*60)) + " seconds");
			 writer.newLine();
			 if(n == 0) {
				 writer.write("Average Number of Iterations: " + String.valueOf(iter));
			 } else {
				 writer.write("Average Number of Iterations: " + String.valueOf(iter/n)); 
			 }
			 writer.newLine();
			 writer.write("Energy of minimum = " +  e + ", found at step " + stepLow);
		} catch (IOException e1) {
		}
		{
			try {
				if (writer != null)
					writer.close();
			} catch (IOException e1) {
			}
		}
	}
	
	public void coords(String optimiser, int atoms, int n, int seed, ArrayList<Minimum> lowest, String potential, Rotation rotation) {
		
		if(potential == "LJ") {
			try {
				writer = new BufferedWriter(new FileWriter("COORDS_" + optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".xyz", false));
				for(int j = 0; j < lowest.size(); j++) {
					if(j > 0) {
						writer.newLine();
					}
					writer.write(String.valueOf(atoms));
					writer.newLine();
					writer.write("Energy of minimum = " +  lowest.get(j).getE() + " found at step " + lowest.get(j).getStep());
					writer.newLine();
					for(int k = 0; k < lowest.get(j).getX().length; k+=3) {
						writer.write("C " + String.valueOf(lowest.get(j).getX()[k]) + " " + String.valueOf(lowest.get(j).getX()[k+1]) + " " + String.valueOf(lowest.get(j).getX()[k+2]) + " ");
						writer.newLine();
					}
				}
			} catch (IOException e1) {
			}
			{
				try {
					if (writer != null)
						writer.close();
				} catch (IOException e1) {
				}
			}
		} else {
			TIP4P tip = new TIP4P();
			double[][] rbcoords = new double[atoms*3][3];
			try {
				writer = new BufferedWriter(new FileWriter("COORDS_" + optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".xyz", false));
				for(int j = 0; j < lowest.size(); j++) {
					if(j > 0) {
						writer.newLine();
					}
					writer.write(String.valueOf(atoms*3));
					writer.newLine();
					writer.write("Energy of minimum = " +  lowest.get(j).getE() + " found at step " + lowest.get(j).getStep());
					writer.newLine();
					rbcoords = tip.viewTIP(lowest.get(j).getX(), rotation);
					
					for(int a = 0; a < rbcoords.length; a++) {
						if(a == 0 || a%3 == 0) {
							writer.write("O " + String.valueOf(rbcoords[a][0]) + " " + String.valueOf(rbcoords[a][1]) + " " + String.valueOf(rbcoords[a][2]) + " ");
							writer.newLine();
						} else {
							writer.write("H " + String.valueOf(rbcoords[a][0]) + " " + String.valueOf(rbcoords[a][1]) + " " + String.valueOf(rbcoords[a][2]) + " ");
							writer.newLine();
						}
					}
				}
			} catch (IOException e1) {
			}
			{
				try {
					if (writer != null)
						writer.close();
				} catch (IOException e1) {
				}
			}
			
		}
	}
	
	public void coordsFilm(String optimiser, int atoms, int n, int seed, double[] x, String potential, Rotation rotation) {
		
		if(potential == "LJ") {
			try {
				writer = new BufferedWriter(new FileWriter("COORDS_FILM_" + optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".xyz", true));
				writer.write(String.valueOf(atoms));
				writer.newLine();
				for(int i = 0; i < x.length; i+=3) {
					writer.write("C " + String.valueOf(x[i]) + " " + String.valueOf(x[i+1]) + " " + String.valueOf(x[i+2]) + " ");
					writer.newLine();
				}
			} catch (IOException e1) {
			}
			{
				try {
					if (writer != null)
						writer.close();
				} catch (IOException e1) {
				}
			}
		} else {
			TIP4P tip = new TIP4P();
			double[][] rbcoords = new double[atoms*3][3];
			try {
				writer = new BufferedWriter(new FileWriter("COORDS_" + optimiser + "_" + atoms + "_" + n + "Steps_" + seed + ".xyz", true));
				writer.write(String.valueOf(atoms*3));
				writer.newLine();
				rbcoords = tip.viewTIP(x, rotation);
				for(int a = 0; a < rbcoords.length; a++) {
					if(a == 0 || a%3 == 0) {
						writer.write("O " + String.valueOf(rbcoords[a][0]) + " " + String.valueOf(rbcoords[a][1]) + " " + String.valueOf(rbcoords[a][2]) + " ");
						writer.newLine();
					} else {
						writer.write("H " + String.valueOf(rbcoords[a][0]) + " " + String.valueOf(rbcoords[a][1]) + " " + String.valueOf(rbcoords[a][2]) + " ");
						writer.newLine();
					}
				}
			} catch (IOException e1) {
			}
			{
				try {
					if (writer != null)
						writer.close();
				} catch (IOException e1) {
				}
			}
			
		}
	}
	
	public void sFilm(int atoms, double[] x, String potential, Rotation rotation, int iter) {
		
		if(potential == "LJ") {
			try {
				writer = new BufferedWriter(new FileWriter("COORDS_FILM_LBFGS" + "_" + atoms + "_LJ_atoms" + ".xyz", true));
				writer.write(String.valueOf(atoms));
				writer.newLine();
				writer.write("Iteration: " + String.valueOf(iter));
				writer.newLine();
				for(int i = 0; i < x.length; i+=3) {
					writer.write("F " + String.valueOf(x[i]) + " " + String.valueOf(x[i+1]) + " " + String.valueOf(x[i+2]) + " ");
					writer.newLine();
				}
			} catch (IOException e1) {
			}
			{
				try {
					if (writer != null)
						writer.close();
				} catch (IOException e1) {
				}
			}
		}
	}
	
	public void fFilm(int atoms, double[] x, String potential, Rotation rotation, int iter) {
		
		if(potential == "LJ") {
			try {
				writer = new BufferedWriter(new FileWriter("COORDS_FILM_LBFGS" + "_" + atoms + "_LJ_atoms" + ".xyz", true));
				writer.write(String.valueOf(atoms));
				writer.newLine();
				writer.write("Iteration: " + String.valueOf(iter));
				writer.newLine();
				for(int i = 0; i < x.length; i+=3) {
					writer.write("Br " + String.valueOf(x[i]) + " " + String.valueOf(x[i+1]) + " " + String.valueOf(x[i+2]) + " ");
					writer.newLine();
				}
			} catch (IOException e1) {
			}
			{
				try {
					if (writer != null)
						writer.close();
				} catch (IOException e1) {
				}
			}
		}
	}

}
