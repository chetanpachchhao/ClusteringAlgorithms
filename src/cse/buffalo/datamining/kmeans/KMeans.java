package cse.buffalo.datamining.kmeans;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import cse.buffalo.datamining.hierarchial.DimensionMatrix;
import cse.buffalo.datamining.pca.Pca;

public class KMeans {
	static ArrayList<ArrayList<Double>> changingkList = new ArrayList<ArrayList<Double>>();
	static ArrayList<ArrayList<Double>> originalKList = new ArrayList<ArrayList<Double>>();
	static int globalIteratorCount = 0; 
	static HashMap<Integer, ArrayList<Integer>> clusterWiseGeneId = new HashMap<Integer, ArrayList<Integer>>();
	static DimensionMatrix dm = new DimensionMatrix("dataset2.txt");
	
	public static void kMean(ArrayList<RowStructure> rowStructureArray){
		changingkList = InputK.initializaK(rowStructureArray);
		int[][] PgroundTruth = new int[rowStructureArray.size()][rowStructureArray.size()];
		instantiatePgroundTruth(PgroundTruth);
		
		do{
			globalIteratorCount++; 
			System.out.println("Iteration "+globalIteratorCount);
			//clear Data Structures
			originalKList.clear();
			clusterWiseGeneId.clear();
			//*****************************************************************************8
			originalKList.addAll(changingkList);
			for(RowStructure row : rowStructureArray){
				ArrayList<Double> differenceArray = new ArrayList<Double>(); 
				for(ArrayList<Double> subKList : changingkList){
					differenceArray.add(findDifference(row.geneDimensions, subKList));
				}
				int indexClusterInKList = differenceArray.indexOf(Collections.min(differenceArray)); 
				//update clusterList
				updateClusterWiseGeneId(indexClusterInKList, row.geneId);
			}
			changeKList(rowStructureArray);
			//System.out.println(changingkList);
			//for( ArrayList<Double> temp : changingkList ){
			//	System.out.println(temp);
			//}
			//System.out.println(clusterWiseGeneId);
		}while(!changingkList.equals(originalKList));
		System.out.println("Converged");
		
		for( ArrayList<Double> temp : changingkList ){
			System.out.println(temp);
		}
		
		for(int index : clusterWiseGeneId.keySet())
			System.out.println(clusterWiseGeneId.get(index));
		
		// Construct PgroundTruth Matrix
		constructPgroundTruth(clusterWiseGeneId, PgroundTruth);
		
		// Calculate Jaccard Coefficient
		calculateJaccardCoefficient(dm.getCgroundTruth(), PgroundTruth);
		
		// Calculate Correlation Coefficient
		dm.constructDistanceMatrix();
		Double[][] distanceMatrix = dm.getDistanceMatrix();
		calculateCorrelation(distanceMatrix, PgroundTruth);
		
		// Get Pca ClusterMatrix
		Pca pca = new Pca(clusterWiseGeneId); 
		pca.printOutput();
	}
	
	
	private static void calculateCorrelation(Double[][] distanceMatrix, int[][] PgroundTruth) {
		Double dBar = calculateMean(distanceMatrix);
		Double cBar = calculateMean(PgroundTruth);
		Double numerator = 0.0, denomenator = 0.0;
		Double dPart = 0.0, cPart = 0.0;
		
		// Calculate numerator & denomenator
		for (int i= 0; i < distanceMatrix.length; i++){
			for (int j = 0; j < distanceMatrix.length; j++){
				if (i != j){
					numerator += (distanceMatrix[i][j] - dBar) * (PgroundTruth[i][j] - cBar);
					dPart += Math.pow((distanceMatrix[i][j] - dBar), 2);
					cPart += Math.pow((PgroundTruth[i][j] - cBar), 2);
				}
			}
		}
//		System.out.println("DPart : " + dPart + " CPart : " + cPart);
		denomenator = Math.sqrt(dPart) * Math.sqrt(cPart);
		
		System.out.println("Correlation is : " + numerator / denomenator);
	}

	private static Double calculateMean(int[][] matrix) {
		Double mean = 0.0, sum = 0.0;
		for (int i = 0; i < matrix.length; i++){
			for ( int j = 0; j < matrix.length; j++){
				if(i != j)
					sum += matrix[i][j];
			}
		}
		mean = (sum * 1.0) / (matrix.length * matrix.length);
		return mean;
	}

	private static Double calculateMean(Double[][] matrix) {
		Double mean = 0.0, sum = 0.0;
		for (int i = 0; i < matrix.length; i++){
			for ( int j = 0; j < matrix.length; j++){
				if(i != j)
					sum += (Double)matrix[i][j];
			}
		}
		mean = sum / (matrix.length * matrix.length);
		return mean;
	}

	
	private static void instantiatePgroundTruth(int[][] PgroundTruth) {
		for (int i = 0; i < PgroundTruth.length; i++){
			for (int j = 0; j < PgroundTruth[0].length; j++)
				PgroundTruth[i][j] = 0;
		}
	}

	
	private static void calculateJaccardCoefficient(int[][] CgroundTruth,
			int[][] PgroundTruth) {

		int m00 = 0, m01 = 0, m10 = 0, m11 = 0;
		for(int i = 0; i < CgroundTruth.length; i++){
			for(int j = 0; j < CgroundTruth.length; j++){
				if (CgroundTruth[i][j] == 0 && PgroundTruth[i][j] == 0)
					m00++;
				if (CgroundTruth[i][j] == 0 && PgroundTruth[i][j] == 1)
					m01++;
				if (CgroundTruth[i][j] == 1 && PgroundTruth[i][j] == 0)
					m10++;
				if (CgroundTruth[i][j] == 1 && PgroundTruth[i][j] == 1)
					m11++;
			}
		}
//		System.out.println("m00 : " + m00 + " m01 : " + m01 + " m10 : " + m10 + " m11 : " + m11);
		System.out.println("Jaccard Coefficient : " + ((m11 * 1.0) / (m11 + m10 + m01)) );
	
	}

	private static void constructPgroundTruth(
			HashMap<Integer, ArrayList<Integer>> clusterWiseGeneId,
			int[][] PgroundTruth) {

		for(int key : clusterWiseGeneId.keySet()){
			if (clusterWiseGeneId.get(key) != null){
				ArrayList<Integer> geneIdList = clusterWiseGeneId.get(key);
				for(int i : geneIdList){
					for(int j : geneIdList){
						PgroundTruth[i-1][j-1] = 1;
					}
				}
			}
		}
		
//		for (int i = 0; i < PgroundTruth.length; i++){
//			System.out.println(Arrays.toString(PgroundTruth[i]));
//		}
	
	}

	private static void updateClusterWiseGeneId(int indexClusterInKList, int geneId) {
		if(clusterWiseGeneId.containsKey(indexClusterInKList)){
			ArrayList<Integer> tempList = clusterWiseGeneId.get(indexClusterInKList);
			tempList.add(geneId); 
			clusterWiseGeneId.put(indexClusterInKList, tempList); 
		}
		else{
			ArrayList<Integer> tempList = new ArrayList<Integer>(); 
			tempList.add(geneId); 
			clusterWiseGeneId.put(indexClusterInKList, tempList); 
		}
	}
	
	private static void changeKList(ArrayList<RowStructure> rowStructureArray) {
		
		for(int indexClusterInKList : clusterWiseGeneId.keySet()){
			ArrayList<Integer> clusterIds = clusterWiseGeneId.get(indexClusterInKList);
			ArrayList<Double> adderList = new ArrayList<Double>(); 
			for(Integer clusterId : clusterIds){
				ArrayList<Double> geneDimensionsForClusterId = new ArrayList<Double>(); 
				geneDimensionsForClusterId = rowStructureArray.get(clusterId - 1).geneDimensions; 
				for(int count = 0; count < geneDimensionsForClusterId.size(); count++){
					if(adderList.size() <= count){
						adderList.add(geneDimensionsForClusterId.get(count)); 
					}
					else{
						adderList.set(count, adderList.get(count) + geneDimensionsForClusterId.get(count)); 
					}
				}
			}
			//divide by size
			int size = clusterIds.size(); 
			for(int  i =0 ; i < adderList.size(); i++){
				adderList.set(i, adderList.get(i) / size );
			}
			//change k List
			changingkList.set(indexClusterInKList, adderList ); 
		}
	}

	private static Double findDifference(ArrayList<Double> geneDimensions,
			ArrayList<Double> subKList) {
		Double distance = 0.0; 
		for(int i = 0; i < geneDimensions.size(); i++){
			Double d1 = geneDimensions.get(i);
			Double d2 = subKList.get(i);
			Double differenceSquare = Math.pow( (d1-d2), 2); 
			distance = distance + differenceSquare; 
		}	
		return Math.sqrt(distance);
	}
	
}

