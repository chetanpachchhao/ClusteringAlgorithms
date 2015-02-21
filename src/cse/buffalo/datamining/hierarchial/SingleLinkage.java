package cse.buffalo.datamining.hierarchial;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import cse.buffalo.datamining.pca.Pca;

public class SingleLinkage {

	public static void main(String[] args) {
		DimensionMatrix dm = new DimensionMatrix("dataset2.txt");
		Double[][] distanceMatrix = dm.getDistanceMatrix();
		int[] dmin = dm.getDMin();
		int sequenceNumber = 0;
		HashMap<Integer, ArrayList<Integer>> mergeStatus = new HashMap<Integer, ArrayList<Integer>>();
		int[][] PgroundTruth = new int[distanceMatrix.length][distanceMatrix.length];
		
		// Instantiate the mergeStatus list and PgroundTruth
		instantiateMergeStatus(mergeStatus, distanceMatrix.length);
		instantiatePgroundTruth(PgroundTruth);
		
		int N = distanceMatrix.length;
		double INFINITY = Double.POSITIVE_INFINITY;
		
		// Help : http://bit.ly/1s1lw2V
		for (sequenceNumber = 0; sequenceNumber < N-5; sequenceNumber++){
			// Find closest pairs of the two clusters
			int cluster1 = 0;
			for (int i = 0; i < N; i++){
				if (distanceMatrix[i][dmin[i]] < distanceMatrix[cluster1][dmin[cluster1]]) 
					cluster1 = i;
			}
			int cluster2 = dmin[cluster1];
			
			// overwrite row cluster1 with minimum of entries in row cluster1 and cluster2
			for (int j = 0; j < N; j++){
				if (distanceMatrix[cluster2][j] < distanceMatrix[cluster1][j]){ 
					distanceMatrix[cluster1][j] = distanceMatrix[cluster2][j];
					distanceMatrix[j][cluster1] = distanceMatrix[cluster2][j];
				}
			}
			distanceMatrix[cluster1][cluster1] = INFINITY;
	
			// infinity-out old row cluster2 and column cluster2
			for (int i = 0; i < N; i++){
				distanceMatrix[cluster2][i] = INFINITY; 
				distanceMatrix[i][cluster2] = INFINITY;
			}
			
			// update dmin and replace ones that previous pointed to cluster2 to point to cluster1
			for (int j = 0; j < N; j++) {
				if (dmin[j] == cluster2) 
					dmin[j] = cluster1;
				if (distanceMatrix[cluster1][j] < distanceMatrix[cluster1][dmin[cluster1]]) 
					dmin[cluster1] = j;
			} 
			
			// update mergeStatus
			cluster1++; cluster2++;
			
			System.out.println("-----------------------------------------------------------------");
			System.out.println("Cluster1 : " + mergeStatus.get(cluster1));
			System.out.println("Cluster2 : " + mergeStatus.get(cluster2));
			ArrayList<Integer> temp, tempMax;
			if (mergeStatus.containsKey(cluster1)){
				temp = mergeStatus.get(cluster1);
				tempMax = mergeStatus.get(cluster2);
				temp.addAll(tempMax);
			}else{
				temp = new ArrayList<Integer>();
				temp.add(cluster2);
			}
			mergeStatus.put(cluster1, temp);
			mergeStatus.put(cluster2, null);
			System.out.println("Final    : " + mergeStatus.get(cluster1));
		}

		
		displayMergeStatus(mergeStatus);
		
		// Construct PgroundTruth Matrix
		constructPgroundTruth(mergeStatus, PgroundTruth);
		
		// Calculate Jaccard Coefficient
		calculateJaccardCoefficient(dm.getCgroundTruth(), PgroundTruth);
		
		// Calculate Correlation Coefficient
		dm.constructDistanceMatrix();
		distanceMatrix = dm.getDistanceMatrix();
		calculateCorrelation(distanceMatrix, PgroundTruth);
		
		// Get Pca ClusterMatrix
		Pca pca = new Pca(mergeStatus); 
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

	private static void calculateJaccardCoefficient(int[][] CgroundTruth, int[][] PgroundTruth) {
		int m00 = 0, m01 = 0, m10 = 0, m11 = 0;
		for(int i = 0; i < CgroundTruth.length; i++){
			for(int j = 0; j < CgroundTruth[0].length; j++){
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

	private static void instantiatePgroundTruth(int[][] PgroundTruth) {
		for (int i = 0; i < PgroundTruth.length; i++){
			for (int j = 0; j < PgroundTruth[0].length; j++)
				PgroundTruth[i][j] = 0;
		}
	}

	private static void constructPgroundTruth(HashMap<Integer, ArrayList<Integer>> mergeStatus, int[][] PgroundTruth) {
		for(int key : mergeStatus.keySet()){
			if (mergeStatus.get(key) != null){
				ArrayList<Integer> geneIdList = mergeStatus.get(key);
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

	private static void displayMergeStatus(HashMap<Integer, ArrayList<Integer>> mergeStatus) {
		for(int key : mergeStatus.keySet()){
			if (mergeStatus.get(key) != null)
				System.out.println(key + " : " + mergeStatus.get(key));
		}
	}

	private static void instantiateMergeStatus(HashMap<Integer, ArrayList<Integer>> mergeStatus, int numberOfKeys) {
		for (int i = 1; i <= numberOfKeys; i++){
			ArrayList<Integer> temp = new ArrayList<Integer>();
			temp.add(i);
			mergeStatus.put(i, temp);
		}
	}

	private static void displayMatrix(Double[][] distanceMatrix){
		for(int i = 0; i < distanceMatrix.length; i++){
			System.out.println(Arrays.toString(distanceMatrix[i]));
		}
	}
}
