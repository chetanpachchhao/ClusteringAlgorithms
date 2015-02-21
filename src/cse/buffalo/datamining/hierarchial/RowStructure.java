package cse.buffalo.datamining.hierarchial;

import java.util.ArrayList;

public class RowStructure {
	int geneId; 
	int groundTruth; 
	ArrayList<Double> geneDimensions = new ArrayList<Double>(); 
	
	public void display() {
	  System.out.println("geneId: "+this.geneId+" | groundTruth: "+this.groundTruth+" | geneDimensions: "+this.geneDimensions);
	}
	
	public static Double computeDistance(RowStructure rs1, RowStructure rs2){
		Double distance = (double) 0;
		if (!rs1.equals(rs2)){
			ArrayList<Double> dimension1 = rs1.geneDimensions;
			ArrayList<Double> dimension2 = rs2.geneDimensions;
			for(int i = 0; i < dimension1.size(); i++){
				distance += Math.pow( (dimension1.get(i) - dimension2.get(i)), 2 );
			}
			distance = Math.sqrt(distance);
		}
		return distance;
	}
}
