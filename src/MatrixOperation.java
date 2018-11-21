import java.util.ArrayList;
import java.util.HashMap;

/**
 * This class contains the matrix operation methods for the Reed Muller code.
 *
 * @author 160021429
 */
public class MatrixOperation {
	private static final int ZERO = 0;
	private static final int ONE = 1;
	private static final int TWO = 2;

	/**
	 * The aim of this method is to generate the generator matrix for the Reed Muller code.
	 *
	 * @param k the value k of ReedMuller(k, r)
	 * @param r the value r of ReedMuller(k, r)
	 * @param dimension the dimension of the reed muller code
	 * @param length the length of the reed muller code
	 * @return g the generator matrix
	 */
	static int[][] generateG(int k, int r, int dimension, int length) {
		int[][] g = new int[dimension][length];

		for (int i = ZERO; i < length; i++) {
			g[ZERO][i] = ONE;
		}

		if (r > ZERO) {
			int limit = k - ONE;
			int row = ONE;

			for (int i = ONE; i <= r; i++) {
				ArrayList<ArrayList<Integer>> list = makeCombinition(limit, i);

				for (int j = ZERO; j < list.size(); j++) { //use for loop to iterate the array list of combinations.
					ArrayList<Integer> integers = list.get(j);

					for (int x = ZERO; x < length; x++) { //for loop for S(C)
						boolean checker = true;

						for (int y = ZERO; y < integers.size(); y++) {
							int val = (int) Math.pow(TWO, integers.get(y));

							if ((val & x) == ZERO) {
								checker = false;
								break;
							}
						}

						if (checker) {
							g[row][x] = ONE;
						}
					}
					row += ONE; //change the row
				} //for loop (which is for iterating ArrayList<ArrayList<Integer>> list) ends

			} // the outer for loop ends

			//convert the generator matrix to the normal form
			MatrixOperation.convertToNormalForm(dimension, length, g);

			//convert the generator matrix to the standard form
			MatrixOperation.convertToStandardForm(length, dimension, g);

		} //if statement ends

		return g;
	}

	/**
	 * The aim of this method is to swap 2 columns.
	 * @param i1 index of the first column
	 * @param i2 index of the second column
	 * @param dimension the dimension of the code
	 * @param g the generator matrix
	 */
	public static void swapCol(int i1, int i2, int dimension, int[][] g) {
		for (int i = ZERO; i < dimension; i++) {
			int temp = g[i][i1];
			g[i][i1] = g[i][i2];
			g[i][i2] = temp;
		}
	}

	/**
	 * The aim of this method is to convert the generator matrix to the standard form.
	 * @param length the length of the code
	 * @param dimension the dimension of the code
	 * @param g the generator matrix
	 */
	public static void convertToStandardForm(int length, int dimension, int[][] g) {
		for (int i = ZERO; i < dimension; i++) {
			if (g[i][i] != ONE) {
				int col = i + ONE;

				while (col < length) {
					if (g[i][col] != ONE) {
						col += ONE;
					} else {
						break;
					}
				}

				swapCol(col, i, dimension, g); //swap the columns
			}
		}
	}

	/**
	 * The aim of this method is to generate the syndrome table.
	 */
	static void getSyndrome(int length, int distance, int[][] h, HashMap<Integer, ArrayList<Integer>> syndrome) {
		int numOfCols_H = h[ZERO].length;
		int numOfMaxError = (distance - ONE) / TWO;
		int limit = length - ONE;

		for (int x = ONE; x <= numOfMaxError; x++) {
			//make the combinations to generate error bits
			ArrayList<ArrayList<Integer>> combinations = makeCombinition(length - ONE, x);

			for (int i = ZERO; i < combinations.size(); i++) {
				int eVector = 0;
				int[] e = new int[length];

				ArrayList<Integer> c = combinations.get(i);

				for (int j = ZERO; j < c.size(); j++) {
					int err = ONE << (limit - c.get(j));
					eVector += err;
					e[c.get(j)] = ONE;
				}

				int sVector = ZERO;

				for (int y = ZERO; y < numOfCols_H; y++) {
					int xor = ZERO;

					for (int z = ZERO; z < h.length; z++) {
						xor ^= (h[z][y] & e[z]); //calculate the syndrome by multiplying error vector and parity check matrix
					}

					sVector <<= 1;
					sVector += xor;
				}

				ArrayList<Integer> list;

				//check if the syndrome table already has same integer array as a key
				if (syndrome.containsKey(sVector)) {
					list = syndrome.get(sVector);
				} else {
					list = new ArrayList<>();
				}

				list.add(eVector);
				syndrome.put(sVector, list); //add a new row to the syndrome table
			}

		}

		// add the all zero row to the syndrome table
		ArrayList<Integer> list = new ArrayList<>();
		list.add(ZERO);

		syndrome.put(ZERO, list);

	}

	/**
	 * The aim of this method is to swap 2 rows.
	 * @param i1 index of the first row
	 * @param i2 index of the second row
	 * @param length the length of the code
	 * @param g the generator matrix
	 */
	public static void swapRow(int i1, int i2, int length, int[][] g) {
		for (int i = ZERO; i < length; i++) {
			int temp = g[i1][i];
			g[i1][i] = g[i2][i];
			g[i2][i] = temp;
		}
	}

	/**
	 * This method converts the generator matrix to the normal form.
	 * @param dimension the dimension of the code
	 * @param length the length of the code
	 * @param g the generator matrix
	 */
	public static void convertToNormalForm(int dimension, int length, int[][] g) {
		int row = ZERO;
		for (int i = ZERO; i < length; i++) {

			for (int j = row; j < dimension; j++) {
				if (g[j][i] != ZERO) {
					if (row != j) {
						swapRow(row, j, length, g); //swap 2 rows to make the normal form
						j = row;
					}

					for (int x = ZERO; x < dimension; x++) {
						if (g[x][i] != ZERO && x != j) {

							/*
							 * Add the current row to all rows above and below which have 1 in the
							 * same column to make a normal form.
							 */
							for (int y = ZERO; y < length; y++) {
								g[x][y] = g[x][y] ^ g[j][y];
							}

						}
					}

					row++;
					break;
				}
			}

		}
	}

	/**
	 * The aim of this method is to make the combinations for the Reed Muller code.
	 * @param n for n_C_k
	 * @param k for n_C_k
	 * @return the array list that contains the all combinations.
	 */
	static ArrayList<ArrayList<Integer>> makeCombinition(int n, int k) {
		ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();

		if (n < ZERO || n < k) { //check if the n is in the correct range
			return result; //if not, return the empty list
		}

		ArrayList<Integer> item = new ArrayList<Integer>();
		findCombinations(n, k, ZERO, item, result);

		return result;
	}

	/**
	 * Find all combinations and add them in the array list.
	 * @param n for n_C_k
	 * @param k for n_C_k
	 * @param start the starting number
	 * @param item the array list for the combination
	 * @param res the array list that contains all combinations (the result of this method)
	 */
	static void findCombinations(int n, int k, int start, ArrayList<Integer> item, ArrayList<ArrayList<Integer>> res) {
		if (item.size() == k) {
			res.add(new ArrayList<Integer>(item));
			return;
		}

		for (int i = start; i <= n; i++) {
			item.add(i);
			findCombinations(n, k, i + ONE, item, res);
			item.remove(item.size() - ONE);
		}
	}

}
