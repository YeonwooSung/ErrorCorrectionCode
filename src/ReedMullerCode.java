import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

/**
 * The class for the Reed Muller code.
 * @author 160021429
 */
public class ReedMullerCode implements IECC {
	private final int ZERO = 0;
	private final int ONE = 1;
	private final int TWO = 2;
	private boolean invalid;
	private int length;
	private int distance;
	private int dimension;
	private int[][] g; // generator matrix
	private int[][] h; // parity check matrix
	private HashMap<Integer, ArrayList<Integer>> syndrome = new HashMap<>();

	/**
	 * The aim of this constructor is to precompute the generator matrix, parity check matrix, and the syndrome table.
	 * @param k the value k for RM(k, r)
	 * @param r the value r for RM(k, r)
	 */
	public ReedMullerCode(int k, int r) {
		if (k >= r && r >= ZERO) {
			this.length = (int) Math.pow(TWO, k);         //length of the code = 2^k
			this.distance = (int) Math.pow(TWO, (k - r)); //minimum distance = 2^(k-r)
			this.dimension = calculateDimension(k, r);

			g = MatrixOperation.generateG(k, r, this.dimension, this.length); // generate the generator matrix

			int parityMatrixLength = this.length - this.dimension; //calculate the number of columns in the parity check matrix

			// generate the parity check matrix H
			h = new int[this.length][parityMatrixLength];

			// set the correct values to the parity check matrix

			for (int i = ZERO; i < this.dimension; i++) {
				for (int j = ZERO; j < parityMatrixLength; j++) {
					h[i][j] = g[i][this.dimension + j];
				}
			}

			// append the identity matrix at the bottom of the parity check matrix
			for (int i = ZERO; i < parityMatrixLength; i++) {
				h[i + this.dimension][i] = ONE;
			}

			MatrixOperation.getSyndrome(this.length, this.distance, h, syndrome); // generate the syndrome table

		} else {
			System.out.println("Wrong parameter: 0 <= r <= k");
			this.invalid = true;
		}
	}

	/**
	 * This method calculates the dimension of the current reed muller code instance.
	 *
	 * @param k the number k of the reed muller code RM(k, r)
	 * @param r the number r of the reed muller code RM(k, r)
	 * @return the dimension of the current reed muller code instance.
	 */
	private int calculateDimension(int k, int r) {
		int sum = ONE;

		for (int i = ONE; i <= r; i++) {
			if (i != ONE) {
				int limit = k - i;
				int temp = k;

				for (int j = k - 1; j > limit; j--) {
					temp *= j;
				}

				for (int j = TWO; j <= i; j++) {
					temp /= j;
				}
				sum += temp;
			} else {
				sum += k;
			}
		}
		return sum;
	}

	/**
	 * The getter of the length attribute.
	 *
	 * @return the length of the code, the number of bits in each encoded block
	 */
	@Override
	public int getLength() {
		return this.length;
	}

	/**
	 * The getter of the attribute dimension.
	 *
	 * @return the dimension of the code, the number of bits in each plain-text block
	 */
	@Override
	public int getDimension() {
		return this.dimension;
	}

	/**
	 * This method encodes the given bit set with the generator matrix, which is precomputed.
	 *
	 * @param plaintext the plain text that should be encoded
	 * @param len the length of the plain text
	 * @return encoded bit set
	 */
	@Override
	public BitSet encode(BitSet plaintext, int len) {
		if (this.invalid) {
			System.out.println("Cannot encode with this instance!");
			return plaintext;
		}

		int numOfBlocks = (len > this.dimension) ? (
				(len % this.dimension == ZERO) ? len / this.dimension : len / this.dimension + ONE
			) : ONE;

		int totalLength = this.length * numOfBlocks;

		BitSet encoded = new BitSet(totalLength);

		int nonZero = plaintext.nextSetBit(ZERO);
		int blockCount = ZERO;
		int index = ZERO;

		for (int i = 0; i < numOfBlocks; i++) {
			int[] text = new int[this.dimension]; //an integer array to store bits of the plain text in it

			for (int j = 0; j < this.dimension; j++) {
				if ((j + index) == nonZero) { //check if the current index is a non zero bit.
					text[j] = ONE;
					nonZero = plaintext.nextSetBit(nonZero + 1);
				}
			}

			int sum = 0;

			for (int x = 0; x < this.length; x++) {
				for (int y = 0; y < this.dimension; y++) {
					int val = g[y][x] * text[y]; // multiply the plain text and generator matrix
					sum = sum ^ val;
				}

				if (sum != ZERO) {
					encoded.set(x + blockCount);
				}
				sum = ZERO; //reset the value of sum to 0
			}
			blockCount += this.length;
			index += this.dimension;
		}

		return encoded;
	}

	/**
	 * This method decodes the given codes by using the precomputed parity check matrix and syndrome table.
	 * @param codetext the code that should be decoded
	 * @param len the length of the code
	 * @param checkIfUnique True if the decode() should check if the error correction could be done uniquely.
	 * @return decoded bit set
	 */
	private BitSet decode(BitSet codetext, int len, boolean checkIfUnique) {
		if (this.invalid) {
			System.out.println("Cannot decode with this instance!");
			return codetext;
		}

		int numOfBlocks = (len > this.length) ? (
				(len % this.length != 0) ? (len / this.length) + 1 : (len / this.length)
			) : 1;

		int totalLength = this.dimension * numOfBlocks;

		BitSet decoded = new BitSet(totalLength);

		int index = ZERO;
		int nonZero = codetext.nextSetBit(ZERO);
		int blockCount = ZERO;

		for (int i = 0; i < numOfBlocks; i++) {
			int[] code = new int[this.length];

			for (int j = 0; j < this.length; j++) {
				if ((j + blockCount) == nonZero) { //check if the current index is a non zero bit.
					code[j] = ONE;
					nonZero = codetext.nextSetBit(nonZero + 1);
				}
			}

			int sVector = 0;

			// multiply H and code to get the syndrome
			for (int x = ZERO; x < h[ZERO].length; x++) {
				int xor = ZERO;
				for (int y = ZERO; y < this.length; y++) {
					xor ^= (code[y] * h[y][x]);
				}

				sVector <<= 1;
				sVector += xor;
			}

			ArrayList<Integer> syndromes = syndrome.get(sVector);
			if (checkIfUnique) {
				int numOfMaxError = (distance - ONE) / TWO;

				if (numOfMaxError == ZERO) {
					return null; //if the current code instance could not correct any error, return null
				}

				if (syndromes.size() > ONE) {
					return null; //if there are more than one error vectors that has same syndrome value, return null
				}
			}

			try {
				//correct the errors to get the closest code
				int e = syndromes.get(ZERO);
				int limit = this.length - 1;

				for (int j = ZERO; j < this.length; j++) {
					int val = ONE << (limit - j);

					if ((val & e) != ZERO) {
						code[j] ^= ONE;
					}
				}
			} catch (NullPointerException e) {
				//If the RM code instance could correct 0 error, NullPointerException might be occurred.
				e.getMessage();
			}

			for (int j = ZERO; j < this.dimension; j++) {
				if (code[j] != ZERO) {
					decoded.set(j + index);
				}
			}

			blockCount += this.length;
			index += this.dimension;
		}

		return decoded;
	}

	/**
	 * Decodes a vector of coded text of any length, padding it to whole number of blocks with 0 bits
     * and then replacing each block with the plaintext corresponding to a closest codeword (in Hamming distance).
     *
     * @param codetext the binary input
     * @param len the length of the codetext
     * @return the decoded version of plaintext (padded to a whole number of blocks)
	 */
	@Override
	public BitSet decodeAlways(BitSet codetext, int len) {
		return decode(codetext, len, false); //decode the given code
	}

	/**
     * Decodes a vector of coded text of any length, padding it to whole number of blocks with 0 bits
     * and then replacing each block with the plaintext corresponding to the uniquie closest codeword (in Hamming distance)
     * if there isn't one then it throws an exception.
     *
     * @param codetext the binary input
     * @param len the length of the codetext
     * @return the decoded version of plaintext (padded to a whole number of blocks)
     * @throws UncorrectableErrorException if there is no uniquely best decoding
     */
	@Override
	public BitSet decodeIfUnique(BitSet codetext, int len) throws UncorrectableErrorException {
		if (this.invalid) {
			System.out.println("Cannot decode with this instance!");
			return codetext;
		}

		BitSet decoded = decode(codetext, len, true); //decode the given code

		if (decoded != null) {
			return decoded;
		} else {
			throw new UncorrectableErrorException();
		}
	}

	/**
	 * This method tells the user about the information of the ReedMullerCode type object.
	 * Or it will tell the user that the current object is an invalid code object.
	 *
	 * @return a message that identifies the current object
	 */
	public String toString() {
		if (this.invalid) {
			return "Invalid code!";
		}

		return ("<Reed Muller: length(" + this.length + "), dimension(" + this.dimension + ")>");
	}

}
