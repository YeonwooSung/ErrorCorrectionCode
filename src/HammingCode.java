import java.util.ArrayList;
import java.util.BitSet;

/**
 * The class for the hamming code.
 * @author 160021429
 */
public class HammingCode implements IECC {
	private final int ZERO = 0;
	private final int ONE = 1;
	private final int TWO = 2;
	private int dimension;
	private int length;
	private boolean invalid;
	private ArrayList<Integer> parityBits = new ArrayList<>();

	/**
	 * This constructor sets the proper values to the fields.
	 * @param r The number of parity check bits.
	 */
	public HammingCode(int r) {
		if (r < 2) {
			this.invalid = true;
		} else {
			length = (int) Math.pow(2, r) - 1;
			dimension = length - r;

			for (int i = ZERO; i < r; i++) {
				parityBits.add(((int) Math.pow(TWO, i)) - ONE);
			}
		}
	}

	/**
	 * The getter for the length.
	 * @return the length of the code
	 */
	@Override
	public int getLength() {
		return this.length;
	}

	/**
	 * The getter for the dimension.
	 * @return the dimension of the code
	 */
	@Override
	public int getDimension() {
		return dimension;
	}

	/**
	 * Set the proper values to the parity bits.
	 * @param bits The array that contains the code bits.
	 * @return code array with proper value of parity bits.
	 */
	private int[] checkAndSetParity(int[] bits) {
		int total = bits.length;
		int size = parityBits.size();
		int[] parity = new int[size];

		for (int i = ZERO; i < size; i++) {
			int parityIndex = parityBits.get(i);
			int parityCount = parityIndex + ONE;

			int sum = ZERO;
			for (int k = parityIndex; k < total; k += parityCount) {

				int count = ZERO;

				while (count < parityCount && k < total) {
					sum += bits[k];
					count += ONE;
					k += ONE;
				}
			}

			// if the result of parity checking is an odd number, set the parity bit as 1
			if (sum % TWO != ZERO) {
				parity[i] = ONE;
			}
		}

		return parity;
	}

	/**
	 * This method encodes the given plain text by using the hamming code.
	 *
	 * @param plaintext The bit set of plain text
	 * @param len The length of the code
	 * @return encoded bit set
	 */
	@Override
	public BitSet encode(BitSet plaintext, int len) {
		if (this.invalid) {
			return null;
		}

		int numOfBlock = (len > this.dimension) ? (
				(len % this.dimension == 0) ? len / this.dimension : (len / this.dimension) + 1
			) : 1;

		int totalLength = this.getLength() * numOfBlock;

		BitSet encoded = new BitSet(totalLength);

		if (plaintext.isEmpty()) { // check if all bits are zero
			return encoded;
		}

		int limit = this.length;
		int start = ZERO;
		int blockIndex = ZERO;

		int nonZero = plaintext.nextSetBit(ZERO);
		int totalIndex = ZERO;

		while (limit <= totalLength) { //iterate the all plain text blocks

			int[] bits = new int[this.length];

			int parityCount = ZERO;
			int count = ZERO;

			int paritySpace = this.parityBits.get(ZERO);
			int index = ZERO;

			while (start < limit) {
				if (count != paritySpace) { //check if the current bit is a parity bit.

					if (nonZero == totalIndex) { //check if the current bit is a set bit

						bits[index] = ONE;
						encoded.set(start);

						// find the next set bit
						nonZero = plaintext.nextSetBit(nonZero + ONE);
					}
					totalIndex += ONE;

				} else {
					parityCount += ONE;
					if (parityCount != parityBits.size()) {
						paritySpace = this.parityBits.get(parityCount);
					}
				}

				start += 1;
				count += 1;
				index += 1;
			}

			int[] parity = checkAndSetParity(bits);

			// set the proper values to the parity bits
			for (int j = ZERO; j < parity.length; j++) {
				if (parity[j] != ZERO) {
					int parityIndex = this.parityBits.get(j) + blockIndex;
					encoded.set(parityIndex);
				}
			}

			limit += this.length;
			blockIndex += this.length;
		}

		return encoded;
	}

	/**
	 * This method corrects the error by parity checking.
	 * @param code The bit set of code that might be corrected
	 * @param totalLength the total length of the code
	 */
	private void errorCorrection(BitSet code, int totalLength) {
		int nonZero = code.nextSetBit(ZERO);
		int start = ZERO;

		while (start < totalLength) {
			int parityCheckBit = ZERO;
			int[] ia = new int[this.length];

			for (int i = ZERO; i < this.length; i++) {
				int index = start + i;
				if (nonZero == index) {
					ia[i] = 1;
					nonZero = code.nextSetBit(nonZero + ONE);
				}
			}

			for (int i = this.parityBits.size() - 1; i >= ZERO; i--) {
				int parityBit = this.parityBits.get(i);
				int pb = parityBit;
				int sum = ZERO;

				while (pb < this.length) {
					int pCount = ZERO;
					while (pCount <= parityBit) {
						sum += ia[pb];
						pb += ONE;
						pCount += ONE;
					}
					pb += (parityBit + 1);
				}

				parityCheckBit = (parityCheckBit << 1) | (sum % 2);
			}

			if (parityCheckBit != ZERO) {
				code.flip(start + parityCheckBit - 1);
			}

			start += this.length;
		}
	}

	/**
	 * Decodes a vector of coded text of any length, padding it to whole number of blocks with 0 bits
     * and then replacing each block with the plaintext corresponding to A closest codeword (in Hamming distance).
     * @param codetext the binary input
     * @param len the length of the codetext
     * @return the decoded version of plaintext (padded to a whole number of blocks)
	 */
	@Override
	public BitSet decodeAlways(BitSet codetext, int len) {
		if (this.invalid) {
			return null;
		}

		errorCorrection(codetext, len); //check and correct the error in the given code text

		int numOfBlocks = (len > this.length) ? (
				(len % this.length != 0) ? (len / this.length) + 1 : (len / this.length)
		) : 1;

		int totalLength = this.dimension * numOfBlocks;

		int startIndex = ZERO;
		int processedLength = ZERO;
		int nonZero = codetext.nextSetBit(ZERO);
		int index = ZERO;

		BitSet decoded = new BitSet(totalLength);

		while (processedLength < totalLength) {
			int parityIndex = ZERO;
			int paritySpace = this.parityBits.get(parityIndex);

			for (int i = 0; i < this.length; i++) {
				if (paritySpace != i) {
					if (nonZero == startIndex) {
						nonZero = codetext.nextSetBit(nonZero + 1);
						decoded.set(index);
						//System.out.println("bit: " + index + " " + startIndex);
					}
					index += ONE;
				} else {

					if (nonZero == startIndex) {
						nonZero = codetext.nextSetBit(nonZero + 1);
					}

					parityIndex += ONE;
					if (parityIndex < this.parityBits.size()) {
						paritySpace = this.parityBits.get(parityIndex);
					}
				}

				startIndex += ONE;
			}

			processedLength += this.dimension;
		}

		return decoded;
	}

	/**
	 * Decodes a vector of coded text of any length, padding it to whole number of blocks with 0 bits
     * and then replacing each block with the plaintext corresponding to A closest codeword (in Hamming distance).
     *
	 * @param codetext The target code text
	 * @param len The length of the given code text
	 * @return the decoded version of plaintext (padded to a whole number of blocks)
	 * @throws UncorrectableErrorException
	 */
	@Override
	public BitSet decodeIfUnique(BitSet codetext, int len) throws UncorrectableErrorException {
		return this.decodeAlways(codetext, len);
	}

	/**
	 * Returns the suitable string that identifies this instance.
	 * @return identifying string
	 */
	public String toString() {
		if (this.invalid) {
			return "Invalid Hamming Code instance";
		}
		return "length(" + this.length + "), dimenstion(" + this.dimension + ")";
	}
}