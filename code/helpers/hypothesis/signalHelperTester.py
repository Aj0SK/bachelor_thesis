from signalHelper import stringAllignment

assert stringAllignment("ACAATG", "ACAATG") == ("ACAATG", "ACAATG")
assert stringAllignment("ACXAATG", "ACAATG") == ("ACXAATG", "AC-AATG")
assert stringAllignment("TTT", "CCC") == ("---TTT", "CCC---")

from signalHelper import countDashes

assert countDashes("AAAC---AC--", 3) == 1
assert countDashes("AAAC---AC--", 2) == 1
assert countDashes("AAAC----A----", 4) == 2

