package gui;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextArea;

import gop.Hop;

import javax.swing.JComboBox;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.awt.event.ActionEvent;

public class View {
	
    public View() {
        JFrame guiFrame = new JFrame();
        
        Hop s = new Hop();
        
        //make sure the program exits when the frame closes
        guiFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        guiFrame.setTitle("GoGUI");
        guiFrame.setSize(750,250);
        guiFrame.setResizable(false);
      
        //This will center the JFrame in the middle of the screen
        guiFrame.setLocationRelativeTo(null);
        
        //The first JPanel contains a JLabel and JCombobox
        final JPanel comboPanel1 = new JPanel();
        final JPanel comboPanel2 = new JPanel();
        final JPanel comboPanel3 = new JPanel();
        
        //Options for the JComboBox 
        String[] optimisers = {"LBFGS", "CG", "SD"};
        String[] lineSearch = {"Wolfe", "WeakWolfe", "BackTrack"};
        String[] potential = {"LJ", "TIP4P"};
        
        JLabel comboOpt = new JLabel("Optimiser:");
        JComboBox<String> optimise = new JComboBox<String>(optimisers);
        JLabel space1 = new JLabel("          ");
        
        comboPanel1.add(comboOpt);
        comboPanel1.add(optimise);
        comboPanel1.add(space1); 
        
        JLabel comboLS = new JLabel("Line Search:");
        JComboBox<String> ls = new JComboBox<String>(lineSearch);
        JLabel space2 = new JLabel("          ");
        
        comboPanel1.add(comboLS);
        comboPanel1.add(ls);
        comboPanel1.add(space2);
        
        JLabel comboPot = new JLabel("System:");
        JComboBox<String> pot = new JComboBox<String>(potential);
        
        comboPanel1.add(comboPot);
        comboPanel1.add(pot);
        
        JLabel nA = new JLabel("System Size");
        JTextArea n = new JTextArea(1,3);
        n.setBorder(BorderFactory.createBevelBorder(2,Color.GRAY,Color.GRAY));
        n.setBackground(Color.LIGHT_GRAY);
        n.setForeground(Color.BLACK);
        n.setVisible(true);
        n.setText("10");
        JLabel space21 = new JLabel("   ");
        
        comboPanel2.add(nA);
        comboPanel2.add(n);
        comboPanel2.add(space21);
        
        JLabel nhop = new JLabel("Number of Hops");
        JTextArea hop = new JTextArea(1,3);
        hop.setBorder(BorderFactory.createBevelBorder(2,Color.GRAY,Color.GRAY));
        hop.setBackground(Color.LIGHT_GRAY);
        hop.setForeground(Color.BLACK);
        hop.setVisible(true);
        hop.setText("2000");
        JLabel space22 = new JLabel("   ");
        
        comboPanel2.add(nhop);
        comboPanel2.add(hop);
        comboPanel2.add(space22);
        
        JLabel ssize = new JLabel("Hop Size");
        JTextArea steps = new JTextArea(1,3);
        steps.setBorder(BorderFactory.createBevelBorder(2,Color.GRAY,Color.GRAY));
        steps.setBackground(Color.LIGHT_GRAY);
        steps.setForeground(Color.BLACK);
        steps.setVisible(true);
        steps.setText("1.0");
        JLabel space23 = new JLabel("   ");
        
        comboPanel2.add(ssize);
        comboPanel2.add(steps);
        comboPanel2.add(space23);
        
        JLabel seed = new JLabel("Seed");
        JTextArea sd = new JTextArea(1,3);
        sd.setBorder(BorderFactory.createBevelBorder(2,Color.GRAY,Color.GRAY));
        sd.setBackground(Color.LIGHT_GRAY);
        sd.setForeground(Color.BLACK);
        sd.setVisible(true);
        sd.setText("2");
        
        comboPanel2.add(seed);
        comboPanel2.add(sd);
        
        JButton beginOpt = new JButton("Begin Optimisation");
        
        beginOpt.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent event) {
            	
            	int nHops = (int)Math.round(Double.valueOf(hop.getText().trim()));
            	int nAtoms = (int)Math.round(Double.valueOf(n.getText().trim()));
            	int seed = (int)Math.round(Double.valueOf(sd.getText().trim()));
            	double step = Double.valueOf(steps.getText().trim());
            	
            	String optim = optimise.getSelectedItem().toString();
            	String line = ls.getSelectedItem().toString();
            	String pote = pot.getSelectedItem().toString();
            	
            	System.out.println("===========================================================");
            	System.out.println("BEGINNING BASIN-HOPPING\n");
            	s.gOptimise(nHops, nAtoms, 2.0, seed, optim, line, pote, "SXNA", false, step, step);
            	System.out.println("\nBASIN-HOPPING COMPLETED");
            	String user = System.getProperty("user.dir");
            	String cmd1 = "java -jar " + user + "/Jmol.jar " + user + "/COORDS_" + optim + "_" + nAtoms + "_" + nHops + "Steps_" + seed + ".xyz";
            	try {
					Runtime.getRuntime().exec(cmd1);
				} catch (IOException e) {
					e.printStackTrace();
				}
            	System.out.println("===========================================================");
            	
            }
        });
        
        //The JFrame uses the BorderLayout layout manager.
        //Put the two JPanels and JButton in different areas.
        guiFrame.add(comboPanel1, BorderLayout.WEST);
        guiFrame.add(comboPanel2, BorderLayout.NORTH);
        guiFrame.add(comboPanel3, BorderLayout.CENTER);
        guiFrame.add(beginOpt, BorderLayout.SOUTH);
        guiFrame.pack();
        
        //make sure the JFrame is visible
        guiFrame.setVisible(true);
    }
    
}