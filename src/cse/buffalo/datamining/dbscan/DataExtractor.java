package cse.buffalo.datamining.dbscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class DataExtractor {

	public static ArrayList<RowStructure> main(String[] args) {
		ArrayList<RowStructure> rowStructureArray = new ArrayList<RowStructure>(); 
		File file = new File("dataset1.txt"); 
		//File file = new File("temp.txt");
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return rowStructureArray;
	}

}
