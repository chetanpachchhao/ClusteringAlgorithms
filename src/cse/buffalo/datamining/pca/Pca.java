package cse.buffalo.datamining.pca;
import java.util.*;

public class Pca {
	HashMap<Integer, ArrayList<Integer>> inputMap = new HashMap<Integer, ArrayList<Integer>>();  
	int[] output = new int[386];
	
	public Pca(HashMap<Integer, ArrayList<Integer>> input){
		this.inputMap = input;
		constructOutput(inputMap);
	}
	
	public void printOutput(){
		for(int i = 0; i < output.length; i++){
			System.out.print(output[i] + "	");
		}
//		System.out.println(Arrays.toString(output));
	}

	private void constructOutput(HashMap<Integer, ArrayList<Integer>> inputMap) {
		int value = 1;
		for(int key : inputMap.keySet()){
			ArrayList<Integer> temp = inputMap.get(key);
			if (temp != null){
				for(int i : temp){
					output[i-1] = value;
				}
				value++;
			}
		}
	}
}
