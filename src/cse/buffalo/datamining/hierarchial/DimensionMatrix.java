package cse.buffalo.datamining.hierarchial;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

public class DimensionMatrix {
	private File file;
	private ArrayList<RowStructure> rowStructureArray = new ArrayList<RowStructure>();
	private Double[][] distanceMatrix;
	private int[] dmin;
	private int[][] CgroundTruth;
	
	public DimensionMatrix(String name){
		this.file = new File(name);
		constructRowStructureArray();
		constructCgroundTruth();
		constructDistanceMatrix();
	}
	
	public ArrayList<RowStructure> getRowStructureArray(){
		return this.rowStructureArray;
	}
	
	public Double[][] getDistanceMatrix(){
		return this.distanceMatrix;
	}
	
	public int[] getDMin(){
		return this.dmin;
	}
	
	public int[][] getCgroundTruth(){
		return this.CgroundTruth;
	}
	
	// Construct the Ground Truth matrix for given data  
	private void constructCgroundTruth(){
		CgroundTruth = new int[rowStructureArray.size()][rowStructureArray.size()];
		for (int i = 0; i < rowStructureArray.size(); i++){
			for(int j = 0; j < rowStructureArray.size(); j++){
				RowStructure rs1 = rowStructureArray.get(i);
				RowStructure rs2 = rowStructureArray.get(j);
				if (rs1.groundTruth == rs2.groundTruth)
					CgroundTruth[i][j] = 1;
				else
					CgroundTruth[i][j] = 0;
			}
		}
		
//		for (int i = 0; i < CgroundTruth.length; i++){
//			System.out.println(Arrays.toString(CgroundTruth[i]));
//		}
	}
	
	public void constructRowStructureArray(){
		try {
			FileReader fr = new FileReader(file); 
			BufferedReader br = new BufferedReader(fr);
			String line = null;
			if(file.exists())
			{
				while ((line = br.readLine()) != null) {
					String[] inputArray = line.split("\\s+");
					RowStructure rowStructure = new RowStructure(); 
					rowStructure.geneId = Integer.parseInt(inputArray[0]); 
					rowStructure.groundTruth = Integer.parseInt(inputArray[1]); 
					for(int i = 2; i < inputArray.length; i++){
						rowStructure.geneDimensions.add(Double.parseDouble(inputArray[i]));
					}
					rowStructureArray.add(rowStructure); 
				}
				fr.close();
				br.close();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void constructDistanceMatrix(){
		double INFINITY = Double.POSITIVE_INFINITY;
		int N = rowStructureArray.size();
		distanceMatrix = new Double[N][N];
		dmin = new int[N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (i == j){ 
					distanceMatrix[i][j] = INFINITY;
				}
				else {
					RowStructure rs1 = rowStructureArray.get(i);
					RowStructure rs2 = rowStructureArray.get(j);
					distanceMatrix[i][j] = RowStructure.computeDistance(rs1, rs2);
				}
				if (distanceMatrix[i][j] < distanceMatrix[i][dmin[i]]) 
					dmin[i] = j;
			}
		}
	}
	
}
