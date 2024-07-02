import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

public class ColorCodingMain {

    /**maxConv Methode aus dem Cygan Paper Sektion 6.2
     * @param a erstes array
     * @param b zweites array
     * @param length länge des ergebnisarrays
     * @return convolution der beiden arrays a und b
     */
    private static int[] maxConvDomain(int[] a, int[] b, int length) {
		int[] convResult = new int[length];
		for (int i = 0; i < length; i++) {
			int max = 0;
			for (int j = 0; j <= i; j++) {
				// c(i) = max i=j {a(j)+b(i-j)}
				try {
					int max2 = a[j] + b[i - j];
					if (max < max2) {
						max = max2;
					}
				} catch (Exception e) {
					convResult[i] = max;
					return convResult;
				}
			}
			convResult[i] = max;
		}

		return convResult;
	}

	public static int[] maxConv(int[] a, int[] b) {
        int n = a.length;

        // Special handling for empty arrays
        if (n == 0) {
            return new int[0];
        }

        int[] result = new int[n];

        for (int k = 0; k < n; k++) {
            for (int i = Math.max(0, k - n + 1); i <= Math.min(k, n - 1); i++) {
                int j = k - i; 
                result[k] = Math.max(result[k], a[i] + b[j]);
            }
        }

        return result;
    }

	private static Item[] shuffleArray(Item[] Z){
		ArrayList<Item> ZList = new ArrayList<Item>();
		Item[] shuffledZ = new Item[Z.length];
		for(int i = 0; i<Z.length;i++){
			ZList.add(Z[i]);
		}
		Collections.shuffle(ZList);

		for(int i = 0; i< ZList.size();i++){
			shuffledZ[i] = ZList.get(i);
		}

		return shuffledZ;
	}

    public static int[][] randomPartition(Item[] Z, int k, int t) {
		int[][] result = new int[k][t]; // z1 bis zk²
		Random r = new Random();
		
		for(int i = 0; i<Z.length;i++){
			result[r.nextInt(k)][Z[i].getWeight()] = Z[i].getValue();
		}
        return result;
    }

	private static int[] maxConvolutionBinaryTree(int[][] array, int t) { // code von ChatGPT
		if (array == null || array.length == 0) {
			return new int[0];
		}

		// Convert the 2D array into a List of 1D arrays
		ArrayList<int[]> arrays = new ArrayList<>();
		for (int[] row : array) {
			arrays.add(row);
		}

		// Perform binary tree convolution until there's only one array left
		while (arrays.size() > 1) {
			ArrayList<int[]> newArray = new ArrayList<>();
			for (int i = 0; i < arrays.size(); i += 2) {
				if (i + 1 < arrays.size()) {
					int[] a = arrays.get(i);
					int[] b = arrays.get(i + 1);
					int[] result = maxConvDomain(a, b, t);
					newArray.add(result);
				} else {
					// If the last array is odd, add it to the new array list
					newArray.add(arrays.get(i));
				}
			}
			arrays = newArray;
		}

		// At this point, arrays contains only one array which is the result
		return arrays.get(0);
	}

	private static int nextPowerOf2(double n) {
        // Check if the number is positive
        if (n <= 0) {
            throw new IllegalArgumentException("Input must be a positive number.");
        }

        // Convert double to an integer
        int intPart = (int) Math.ceil(n);

        // If intPart is already a power of 2, return intPart
        if ((intPart & (intPart - 1)) == 0) {
            return intPart;
        }

        // Calculate the next power of 2
        intPart--;
        intPart |= intPart >> 1;
        intPart |= intPart >> 2;
        intPart |= intPart >> 4;
        intPart |= intPart >> 8;
        intPart |= intPart >> 16;
        intPart++;

        return intPart;
    }

	private static int[] colorCoding(Item[] Z, int t, int k, double delta){
		if(Z.length == 0){
			return new int[t];
		}
		int log = (int)(Math.ceil((Math.log(1.0/delta))/(Math.log(4.0/3.0))));
		int[][] P = new int[log][];
		for(int i = 1; i<log;i++){
			//randomly partition Z into k^2 pieces where Z[i] is value of Item i
			int[][] partOfZ = new int[k*k][t];
			shuffleArray(Z);
			partOfZ = randomPartition(Z, k*k, t);
			P[i] = maxConvolutionBinaryTree(partOfZ, t);
		}
		int[] W = new int[t];
		for (int j = 1; j<t;j++){
			int max = 0;
			for(int i = 1; i< P.length;i++){
				if(P[i][j] > max){
					max = P[i][j]; // index oder wert??
				}
			}
			W[j] = max;
		}
		return W;
	}
	private static int[] ColorCodingLayer(Item[] L, int t, int i, double delta){
		int l = (int) Math.pow(2, i);
		if(l<(Math.log(l/delta)/Math.log(2))){
			return colorCoding(L,t,l,delta);
		}
		int m = nextPowerOf2(l/(Math.log(l/delta)/Math.log(2)));
		ArrayList<Item>[] A = new ArrayList[m+1];
		for (int a = 0; a <= m; a++) { 
            A[a] = new ArrayList<Item>(); 
        } 
		//ArrayList<ArrayList<Item>> A = new ArrayList<>();
		Random r = new Random();
		for(int p = 0; p<L.length;p++){
			A[1 + r.nextInt(m)].add(L[p]);
		}
		int gamma = (int)(Math.round((6*(Math.log(l/delta)/Math.log(2)))));
		int capacity = 2 * gamma * t / l;
		if (capacity > t) {
			capacity = t;
		}

		int[][] P = new int[m][];
		for(int j = 1;j<m;j++){
			Item[] aArray = A[j].toArray(new Item[A[j].size()]);
			P[j] = colorCoding(aArray, capacity, gamma, delta/l); // sehr hohe zahlen
		}
		for(int h = 1;h<(Math.log(m)/Math.log(2));h++){
			for(int j = 1;j<m/Math.pow(2, h);j++){
				P[j] = maxConvDomain(P[(2*j)-1], P[2*j], (int)Math.pow(2, h)*2*gamma*t/l);
			}
		}
		return P[1];
	}

	private static int[] Knapsack(Item[] Z, int t, double delta){
		int log = (int)(Math.ceil(Math.log(Z.length)/Math.log(2)));
		Item[][] Layers = new Item[log+1][]; // log+1 weil an der logten stelle gibt es auch ein eintrag
		for (int i = 1; i <= log; i++) {
			double under = t/Math.pow(2,i);
			double upper = t/Math.pow(2, i-1);
			ArrayList<Item> layer = new ArrayList<>();
			if (i == log) {
				under = -1; 
				upper = t / Math.pow(2, log-1);
			}
			for(int j=0;j<Z.length;j++){
				if(Z[j].getWeight() > under && Z[j].getWeight()< upper){
					layer.add(Z[j]);
				}
			}

			Item[] arrayLayer = new Item[layer.size()];
			for(int j = 0; j<layer.size();j++){
				arrayLayer[j] = layer.get(j);
			}
			Layers[i] = arrayLayer;
		}
		int[] W = new int[t];
		for(int i = 1;i<=log;i++){
			int[] P_i = ColorCodingLayer(Layers[i], t, i, delta/log);
			W = maxConv(W, P_i);
		}
		return W;
	}

	/* private static HashMap<Integer,Integer> n_and_wmax(String name) throws FileNotFoundException{
		HashMap<Integer,Integer> result = new HashMap<Integer,Integer>();
		Scanner sc = new Scanner(new File(name));
		result.put(Integer.parseInt(sc.next()),Integer.parseInt(sc.next()));
		sc.close();
		return result;
	} */

	private static Item[] readFile(String name) throws IOException {
		Item[] result = null;
		try {
			File myObj = new File(name);
			Scanner sc = new Scanner(myObj);
			int m = Integer.parseInt(sc.next());
			int max = Integer.parseInt(sc.next());
			result = new Item[m];
			for (int i = 0; i < m; i++) {
				int profit = Integer.parseInt(sc.next());
				int weight = Integer.parseInt(sc.next());
				result[i] = new Item(profit, weight);
			}
			sc.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return result;
	}

	private static int readMaxCapacity(String name) throws FileNotFoundException{
		int maxCap = 0;
		Scanner sc = new Scanner(new File(name));
		if (sc.hasNext()){
			sc.next();
			if(sc.hasNext()){
				maxCap = Integer.parseInt(sc.next());
			}
		}
		sc.close();
		return maxCap;
	}

	private static int readOptimumResult(String name) throws IOException {
		File myObj = new File(name);
		Scanner sc = new Scanner(myObj);
		int result = Integer.parseInt(sc.nextLine().trim());
		sc.close();
		return result;
	}

	private static void writeComparisonToLog(int optimumResult, int[] knapsackResult, long duration, String logFile, String setname) {
		try {
			FileWriter writer = new FileWriter(logFile);
			writer.write("Set: " + setname + "\n");
			writer.write("Optimum Result: " + optimumResult + "\n");
			writer.write("Knapsack Result: ");
			for(int i = 0; i < knapsackResult.length; i++){
				writer.write(knapsackResult[i] + ", ");
			}
			
			writer.write("\nDuration: " + duration + " milliseconds\n");
			// Vergleiche die Ergebnisse
			writer.write("Die Ergebnisse liegen " + (optimumResult - knapsackResult[knapsackResult.length-1]) + " auseinander");
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
    public static void main(String[] args) {
		String user = "C:\\Users\\levan\\Desktop\\knapsack\\Knapsack\\";
		String[] itemsLargeScale = { "large_scale\\knapPI_1_100_1000_1",
				"large_scale\\knapPI_1_200_1000_1",
				"large_scale\\knapPI_1_500_1000_1",
				"large_scale\\knapPI_1_1000_1000_1",
				"large_scale\\knapPI_1_2000_1000_1",
				"large_scale\\knapPI_1_5000_1000_1",
				"large_scale\\knapPI_1_10000_1000_1",
				"large_scale\\knapPI_2_100_1000_1",
				"large_scale\\knapPI_2_200_1000_1",
				"large_scale\\knapPI_2_500_1000_1",
				"large_scale\\knapPI_2_1000_1000_1",
				"large_scale\\knapPI_2_2000_1000_1",
				"large_scale\\knapPI_2_5000_1000_1",
				"large_scale\\knapPI_2_10000_1000_1",
				"large_scale\\knapPI_3_100_1000_1",
				"large_scale\\knapPI_3_200_1000_1",
				"large_scale\\knapPI_3_500_1000_1",
				"large_scale\\knapPI_3_1000_1000_1",
				"large_scale\\knapPI_3_2000_1000_1", };

		String[] optimaLargeScale = {
				"large_scale-optimum\\knapPI_1_100_1000_1",
				"large_scale-optimum\\knapPI_1_200_1000_1",
				"large_scale-optimum\\knapPI_1_500_1000_1",
				"large_scale-optimum\\knapPI_1_1000_1000_1",
				"large_scale-optimum\\knapPI_1_2000_1000_1",
				"large_scale-optimum\\knapPI_1_5000_1000_1",
				"large_scale-optimum\\knapPI_1_10000_1000_1",
				"large_scale-optimum\\knapPI_2_100_1000_1",
				"large_scale-optimum\\knapPI_2_200_1000_1",
				"large_scale-optimum\\knapPI_2_500_1000_1",
				"large_scale-optimum\\knapPI_2_1000_1000_1",
				"large_scale-optimum\\knapPI_2_2000_1000_1",
				"large_scale-optimum\\knapPI_2_5000_1000_1",
				"large_scale-optimum\\knapPI_2_10000_1000_1",
				"large_scale-optimum\\knapPI_3_100_1000_1",
				"large_scale-optimum\\knapPI_3_200_1000_1",
				"large_scale-optimum\\knapPI_3_500_1000_1",
				"large_scale-optimum\\knapPI_3_1000_1000_1",
				"large_scale-optimum\\knapPI_3_2000_1000_1", };

		String[] itemsLowScale = { "low-dimensional\\f1_l-d_kp_10_269",
				"low-dimensional\\f2_l-d_kp_20_878",
				"low-dimensional\\f3_l-d_kp_4_20",
				"low-dimensional\\f4_l-d_kp_4_11",
				"low-dimensional\\f5_l-d_kp_15_375",
				"low-dimensional\\f6_l-d_kp_10_60",
				"low-dimensional\\f7_l-d_kp_7_50",
				"low-dimensional\\f8_l-d_kp_23_10000",
				"low-dimensional\\f9_l-d_kp_5_80",
				"low-dimensional\\f10_l-d_kp_20_879" };

		String[] optimaLowScale = {
				"low-dimensional-optimum\\f1_l-d_kp_10_269",
				"low-dimensional-optimum\\f2_l-d_kp_20_878",
				"low-dimensional-optimum\\f3_l-d_kp_4_20",
				"low-dimensional-optimum\\f4_l-d_kp_4_11",
				"low-dimensional-optimum\\f5_l-d_kp_15_375",
				"low-dimensional-optimum\\f6_l-d_kp_10_60",
				"low-dimensional-optimum\\f7_l-d_kp_7_50",
				"low-dimensional-optimum\\f8_l-d_kp_23_10000",
				"low-dimensional-optimum\\f9_l-d_kp_5_80",
				"low-dimensional-optimum\\f10_l-d_kp_20_879" };

		String logFile = "C:\\Users\\levan\\Desktop\\knapsack\\Knapsack\\Log.txt"; // Pfad zur Logdatei

		try {
			String[] testSet = {};
			String[] testOptimum = {};
			System.out.println("1 für Large Scale | 2 für Low Scale: ");
			Scanner scale = new Scanner(System.in);
			int set = scale.nextInt();

			if (set == 1) {
				testSet = itemsLargeScale;
				testOptimum = optimaLargeScale;
				System.out.println(
						"Datensatz wählen: \n1: knapPI_1_100_1000_1 \n2: knapPI_1_200_1000_1 \n3: knapPI_1_500_1000_1 "
								+ "\n4: knapPI_1_1000_1000_1 \n5: knapPI_1_2000_1000_1 \n6: knapPI_1_5000_1000_1 \n7: knapPI_1_10000_1000_1 \n8: knapPI_2_100_1000_1 \n9: knapPI_2_200_1000_1 "
								+ "\n10: knapPI_2_500_1000_1 \n11: knapPI_2_1000_1000_1 \n12: knapPI_2_2000_1000_1 \n13: knapPI_2_5000_1000_1 "
								+ "\n14: knapPI_2_10000_1000_1 \n15: knapPI_3_100_1000_1 \n16: knapPI_3_200_1000_1 \n17: knapPI_3_500_1000_1 \n18: knapPI_3_1000_1000_1 \n19: knapPI_3_2000_1000_1 ");
			}
			if (set == 2) {
				testSet = itemsLowScale;
				testOptimum = optimaLowScale;
				System.out.println("Datensatz wählen: \n1: f1_l-d_kp_10_269 \n2: f2_l-d_kp_20_878 \n3: f3_l-d_kp_4_20 "
						+ "\n4: f4_l-d_kp_4_11 \n5: f5_l-d_kp_15_375 \n6: f6_l-d_kp_10_60 \n7: f7_l-d_kp_7_50 \n8: f8_l-d_kp_23_10000 \n9: f9_l-d_kp_5_80 "
						+ "\n10: f10_l-d_kp_20_879");
			}
			Scanner dataSet = new Scanner(System.in);
			int input = dataSet.nextInt() - 1;
			
			System.out.println("Deltawert zwischen [0,1] wählen: ");
			Scanner sc = new Scanner(System.in);
			double delta = sc.nextDouble();
			scale.close();
			dataSet.close();
			sc.close();


			// Lese die Items aus der Datei ein
			Item[] items = readFile(user + testSet[input]);

			// Lese das Optimumergebnis aus der Datei ein
			int optimumResult = readOptimumResult(user + testOptimum[input]);

			// Starte die Zeitmessung
			long startTime = System.nanoTime();

			// Führe den Knapsack-Algorithmus aus, um das Ergebnis zu erhalten
			int[] knapsackResult = Knapsack(items, readMaxCapacity(user + testSet[input]), delta); // array mit profiten alle werte off siehe paragraph nach knapsack algorithm im paper

			// Beende die Zeitmessung
			long endTime = System.nanoTime();
			long duration = (endTime - startTime) / 1000000; // Zeit in Millisekunden

			// Vergleiche die Ergebnisse und schreibe das Ergebnis in die Logdatei
			writeComparisonToLog(optimumResult, knapsackResult, duration, logFile, testSet[input]);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
