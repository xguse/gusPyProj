from gusPyCode.defs.bioDefs import ParseFastA

fasPath = '/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.masked.fas'

fasParser = ParseFastA(fasPath)

rec     = fasParser.getNext()
recDict = fasParser.toDict()
None