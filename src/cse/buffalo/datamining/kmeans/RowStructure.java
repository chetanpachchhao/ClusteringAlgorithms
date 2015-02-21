package cse.buffalo.datamining.kmeans;

import java.util.ArrayList;

public class RowStructure {
	int geneId; 
	int groundTruth; 
	ArrayList<Double> geneDimensions = new ArrayList<Double>(); 
	public void display() {
		  System.out.println("geneId: "+this.geneId+" | groundTruth: "+this.groundTruth+" | geneDimensions: "+this.geneDimensions);
	 }
}
