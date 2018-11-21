/**
 * The factory that creates the ReedMullerCode instance or the HammingCode instance.
 * @author 160021429
 */
public class ECCFactory implements IECCFactory {

	/**
	 * This method makes the HammingCode instance.
	 * @param r The number of parity check bits.
	 * @return the HammingCode instance
	 */
	@Override
	public IECC makeHammingCode(int r) {
		IECC hammingCode = new HammingCode(r);
		return hammingCode;
	}
	/**
	 * This method makes the ReedMullerCode instance.
	 * @param k for ReedMullerCode(k, r)
	 * @param r for ReedMullerCode(k, r)
	 * @return the ReedMullerCode instance
	 */
	@Override
	public IECC makeReedMullerCode(int k, int r) {
		IECC reedMullerCode = new ReedMullerCode(k, r);
		return reedMullerCode;
	}

}
