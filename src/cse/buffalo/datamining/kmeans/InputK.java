package cse.buffalo.datamining.kmeans;

import java.util.ArrayList;

public class InputK {
	//rowIdsForKCentroids	
	static int[] rowIdsForKCentroids = new int[]{0, 67, 201, 276, 330};
	//static int[] rowIdsForKCentroids = new int[]{0, 3, 6}; 
	public static ArrayList<ArrayList<Double>> initializaK(ArrayList<RowStructure> rowStructureArray){
		ArrayList<ArrayList<Double>> kList = new ArrayList<ArrayList<Double>>(); 
		//For start let us take first 5 geneId as centroids
		for(int i = 0; i < rowIdsForKCentroids.length; i++){
			ArrayList<Double> subKList = new ArrayList<Double>(); 
			for(Double tempVal : rowStructureArray.get(rowIdsForKCentroids[i]).geneDimensions){
				subKList.add(tempVal);
			}
			kList.add(subKList);
		}
		return kList;
	}
}
