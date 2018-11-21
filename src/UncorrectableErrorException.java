/**
 * The custom exception class.
 *
 * @author 160021429
 */
public class UncorrectableErrorException extends Exception {
	/**
	 * The default serial number.
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * The default constructor, which calls the super class's constructor.
	 */
	public UncorrectableErrorException() {
		super();
	}

	/**
	 * The constructor with one argument, which calls the super class's constructor.
	 * @param message The error message
	 */
	public UncorrectableErrorException(String message) {
		super(message);
	}
}
