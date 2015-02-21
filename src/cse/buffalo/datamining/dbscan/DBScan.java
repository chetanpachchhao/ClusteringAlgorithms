package cse.buffalo.datamining.dbscan;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import cse.buffalo.datamining.hierarchial.DimensionMatrix;
import cse.buffalo.datamining.pca.Pca;

public class DBScan {
	
	ArrayList<RowStructure> rowStructureArray = new ArrayList<RowStructure>();
	ArrayList<Integer> visited = new ArrayList<Integer>();
	HashMap<Integer, ArrayList<Integer>> clusters=new HashMap<>();
	public static Double EPS=0.59;//1.95;
	public static Integer minimumPoints=4;//2;
	
	public DBScan(ArrayList<RowStructure> rowStructureArray){
		this.rowStructureArray=rowStructureArray;
		for(int i=0;i<rowStructureArray.size();i++)
			visited.add(0);
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
    DimensionMatrix dm = new DimensionMatrix("dataset1.txt");
    Double[][] distanceMatrix = dm.getDistanceMatrix();
    int[][] PgroundTruth = new int[distanceMatrix.length][distanceMatrix.length];
		
		DBScan db1=new DBScan(DataExtractor.main(null));

		db1.DBScan(db1,EPS,minimumPoints);
		System.out.println("Printing Clusters");
		for(Map.Entry<Integer, ArrayList<Integer>> entry : db1.clusters.entrySet()){
			System.out.println("Cluster-> "+entry+"    size:"+entry.getValue().size());
		}
		
    
    // Construct PgroundTruth Matrix
    constructPgroundTruth(db1.clusters, PgroundTruth);
    
    // Calculate Jaccard Coefficient
    calculateJaccardCoefficient(dm.getCgroundTruth(), PgroundTruth);
    
    // Calculate Correlation Coefficient
    dm.constructDistanceMatrix();
    distanceMatrix = dm.getDistanceMatrix();
    calculateCorrelation(distanceMatrix, PgroundTruth);
    
		// Get Pca ClusterMatrix
		Pca pca = new Pca(db1.clusters); 
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
//    System.out.println("DPart : " + dPart + " CPart : " + cPart);
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
//    System.out.println("m00 : " + m00 + " m01 : " + m01 + " m10 : " + m10 + " m11 : " + m11);
    System.out.println("Jaccard Coefficient : " + ((m11 * 1.0) / (m11 + m10 + m01)) );
  }

  
  private static void constructPgroundTruth(HashMap<Integer, ArrayList<Integer>> mergeStatus, int[][] PgroundTruth) {
    for(int key : mergeStatus.keySet()){
      if (mergeStatus.get(key) != null){
        ArrayList<Integer> geneIdList = mergeStatus.get(key);
        if (key == -1){
          for(int i:geneIdList){
            PgroundTruth[i-1][i-1] = 1;
          }
        }else{
          for(int i : geneIdList){
            for(int j : geneIdList){
              PgroundTruth[i-1][j-1] = 1;
            }
          }
        }
      }
    }
    
//    for (int i = 0; i < PgroundTruth.length; i++){
//      System.out.println(Arrays.toString(PgroundTruth[i]));
//    }
  }


	public void DBScan(DBScan dbscan, Double epsilon,Integer minPts) {
		System.out.println("in dbsacn");
		Integer currentCluster=0;
		int count=0;
		for(int i=0;i<dbscan.rowStructureArray.size();i++){
			ArrayList<Integer> neighborPoints=new ArrayList<>();
			if(dbscan.visited.get(i)==0){
				dbscan.visited.set(i, 1);
				neighborPoints=regionQuery(dbscan.rowStructureArray.get(i), epsilon);
				if(neighborPoints.size()<minPts){
					ArrayList<Integer> pointsInCurrentCluster=new ArrayList<>();
					if(clusters.containsKey(-1)){
						pointsInCurrentCluster=clusters.get(-1);
						pointsInCurrentCluster.add(dbscan.rowStructureArray.get(i).geneId);
						clusters.put(-1, pointsInCurrentCluster);
					}
					else{
						pointsInCurrentCluster=new ArrayList<>();
						pointsInCurrentCluster.add(dbscan.rowStructureArray.get(i).geneId);
						clusters.put(-1, pointsInCurrentCluster);
					}
				}
				else{
					currentCluster++;
					this.expandCluster(dbscan.rowStructureArray.get(i), neighborPoints, currentCluster, epsilon, minPts);
				}
			}	
				
		}
	}
	
	public void expandCluster(RowStructure P,ArrayList<Integer> neighborPoints,Integer currentCluster,Double epsilon,Integer minPts) {
		ArrayList<Integer> pointsInCurrentCluster=new ArrayList<>();
		if(clusters.containsKey(currentCluster)){
			pointsInCurrentCluster=clusters.get(currentCluster);
			pointsInCurrentCluster.add(P.geneId);
			clusters.put(currentCluster, pointsInCurrentCluster);
		}
		else{
			pointsInCurrentCluster.add(P.geneId);
			clusters.put(currentCluster, pointsInCurrentCluster);
		}
		
		
		
		for(int i=0;i<neighborPoints.size();i++){
			for(int j=0;j<this.rowStructureArray.size();j++){
				if(neighborPoints.get(i)==this.rowStructureArray.get(j).geneId){
					if(this.visited.get(j)==0){
						this.visited.set(j, 1);
						ArrayList<Integer> neighborPointsNew=new ArrayList<>();
						neighborPointsNew=regionQuery(this.rowStructureArray.get(j), epsilon);
						if(neighborPointsNew.size()>=minPts){
							for(int k=0;k<neighborPointsNew.size();k++){
								neighborPoints.add(neighborPointsNew.get(k));
							}
						}
					}
					int assignedFlag=0;
					for(int k=1;k<=currentCluster;k++){
						if(clusters.get(k).contains(this.rowStructureArray.get(j).geneId)){
							assignedFlag=1;
							break;
						}
					}
					if(assignedFlag==0){
						pointsInCurrentCluster= new ArrayList<>();
						if(clusters.containsKey(currentCluster)){
							 pointsInCurrentCluster=clusters.get(currentCluster);
							 pointsInCurrentCluster.add(this.rowStructureArray.get(j).geneId);
							 clusters.put(currentCluster, pointsInCurrentCluster);
						}
						else{
							pointsInCurrentCluster.add(this.rowStructureArray.get(j).geneId);
							clusters.put(currentCluster, pointsInCurrentCluster);
						}
					}
				}
			}
		}
	}
	//regionquery
	public ArrayList<Integer> regionQuery(RowStructure P,Double epsilon) {
		int count=0;
		ArrayList<Integer> neighbourPoints=new ArrayList<>();
		for(int i=0;i<this.rowStructureArray.size();i++){
			Double distance=0.0;
				count++;
				//distance = computeDistance(P.geneDimensions, this.rowStructureArray.get(i).geneDimensions);
				distance = computeDistance(P, this.rowStructureArray.get(i));
				if(distance<=EPS)
					neighbourPoints.add(this.rowStructureArray.get(i).geneId);
		}
		return neighbourPoints;
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
