package opt;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import pot.TIP4P;
import rotation.Rotation;

public class FilmWriter {
	
	private BufferedWriter writer = null;
	
	public void coordsFilm(int atoms, double[] x, String potential, Rotation rotation, int iter) {
		
		if(potential == "LJ") {
			try {
				writer = new BufferedWriter(new FileWriter("COORDS_FILM_LBFGS" + "_" + atoms + "_LJ_atoms" + ".xyz", true));
				writer.write(String.valueOf(atoms));
				writer.newLine();
				writer.write("Iteration: " + String.valueOf(iter));
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
				writer = new BufferedWriter(new FileWriter("COORDS_FILM_LBFGS" + "_" + atoms + "_TIP4P_mol" + ".xyz", true));
				writer.write(String.valueOf(atoms*3));
				writer.newLine();
				writer.write("Iteration: " + String.valueOf(iter));
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

}
