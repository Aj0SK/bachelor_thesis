from signalHelper import stringAllignment

assert stringAllignment("ACAATG", "ACAATG") == ("ACAATG", "ACAATG")
assert stringAllignment("ACXAATG", "ACAATG") == ("ACXAATG", "AC-AATG")
assert stringAllignment("TTT", "CCC") == ("---TTT", "CCC---")
