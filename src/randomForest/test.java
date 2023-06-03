package randomForest;

import java.util.Arrays;
import java.util.TreeMap;

public class test {

    public static void main(String[] args) {
        double[] key = {3.3,5.5,2.2,4.2,0.0};
        int[] data = {4,7,5,8,9};
        
        sortArrays(key, data);
        
        for (int i = 0; i < 5; i++) {
            System.out.print(key[i] + ", ");
        }
        System.out.println();
        for (int i = 0; i < 5; i++) {
            System.out.print(data[i] + ", ");
        }

    }
    
    public static void sortArrays(double[] key, int[] data) {
        TreeMap<Double, Integer> helper = new TreeMap<Double, Integer>();
        for (int i = 0; i < key.length; i++) {
            helper.put(key[i], data[i]);
        }
        Arrays.sort(key);
        for (int i = 0; i < data.length; i++) {
            data[i] = helper.get(key[i]);
        }
    }

}
