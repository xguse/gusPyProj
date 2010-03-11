
class InvalidFileFormatException (Exception):
    pass

class CommandLineException (Exception):
    pass

class UnexpectedException(Exception):
    """Used for errors that 'should not happen'.  For example, if you've checked
    that the length is greater than some value and then you check later and the 
    length is shorter, you'd throw this exception."""
    pass

class InvalidInputException(Exception):
    """Raised when a function has been given invalid input.  For example, if a variable
    was expected and None was passed in."""
    pass

class InvalidQuality(Exception):
    """Raised when a function has been an invalid quality value."""
    pass

class InvalidFastq(Exception):
    """Raised when the input fastq format is invalid."""
    pass