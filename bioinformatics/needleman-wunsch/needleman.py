from pprint import pprint

def allMax(*args, key=None):
    """
    Finds all maxes in args.

    :param args:
    :param key: function by which item in arg should be maximum determined
    :return: all maxes from args
    """

    maxes = [args[0], ]

    if key is None:
        for x in args[1:]:
            for y in maxes.copy():
                if x > y:
                    maxes = [x, ]
                    break
                elif x == y:
                    maxes.append(x)
                    break
    else:
        for x in args[1:]:
            z = key(x)
            for y in maxes.copy():
                zz = key(y)
                if z > zz:
                    maxes = [x, ]
                    break
                elif z == zz:
                    maxes.append(x)
                    break
    return maxes

def _operation(s1, s2, al1, al2, x, y, operation):
    if operation == "d":
        al1 += s1[x-1]
        al2 += s2[y-1]
        y -= 1
        x -= 1
        return (al1, al2, x, y)

    if operation == "u":
        al1 += "-"
        al2 += s2[y-1]
        y -= 1
        return (al1, al2, x, y)

    if operation == "l":
        al1 += s1[x-1]
        al2 += "-"
        x -= 1
        return (al1, al2, x, y)

def _backtrackAlts(matrix, s1, s2, x, y, al1="", al2=""):
    """
    Finds all alignments of s1 and s2.

    :param matrix: scoring matrix
    :param s1: first string to align
    :param s2: second string to align
    :param x, y: in first call, the coordinates of right bottom corner of the scoring matrix (i.e. len(s1), len(s2) )
    :return: yields all the possible alignments of s1 and s2 as (al1, al2). First yielded strings are the same as from _backtrack()
    """

    if y == 0 and x == 0:
        yield (al1[::-1], al2[::-1])
    elif len(matrix[y][x]) > 1:
        for alt in matrix[y][x]:
            operation = alt[1]
            _al1, _al2, _x, _y = _operation(s1, s2, al1, al2, x, y, operation)
            yield from _backtrackAlts(matrix, s1, s2, _x, _y, al1=_al1, al2=_al2)
    else:
        operation = matrix[y][x][0][1]
        al1, al2, x, y = _operation(s1, s2, al1, al2, x, y, operation)
        yield from _backtrackAlts(matrix, s1, s2, x, y, al1=al1, al2=al2)

def _backtrack(matrix, s1, s2):
    """
    Finds the first possible alignment of s1 and s2.

    :param matrix: scoring matrix
    :param s1: first string to align
    :param s2: second string to align
    :return: (al1, al2), aligned strings
    """

    al1 = ""
    al2 = ""
    x = len(s1)
    y = len(s2)

    while True:
        if y == 0 and x == 0:
            return (al1[::-1], al2[::-1])
        else:
            operation = matrix[y][x][0][1]
            al1, al2, x, y = _operation(s1, s2, al1, al2, x, y, operation)

def _backtrackRandomized(matrix, s1, s2):
    """
    Finds the random alignments of s1 and s2 if alternative alignments exist.

    :param matrix: scoring matrix
    :param s1: first string to align
    :param s2: second string to align
    :return: (al1, al2), aligned strings
    """

    from random import choice

    al1 = ""
    al2 = ""
    x = len(s1)
    y = len(s2)

    while True:
        if y == 0 and x == 0:
            return (al1[::-1], al2[::-1])
        else:
            operation = matrix[y][x][choice(range(len(matrix[y][x])))][1]
            al1, al2, x, y = _operation(s1, s2, al1, al2, x, y, operation)

def generateMatrix(s1, s2, substMatrix=None, extGapPenalty=False):
    """
    Generates scoring matrix.

    :param s1: first string to align
    :param s2: second string to align
    :param substMatrix: substitution matrix to use. Expected format:
    {
        "A": {"A": +1, "T": -3, "C": -3, "G": -3, "N": -3, },
        "T": {"A": -3, "T": +1, "C": -3, "G": -3, "N": -3, },
        "C": {"A": -3, "T": -3, "C": +1, "G": -3, "N": -3, },
        "G": {"A": -3, "T": -3, "C": -3, "G": +1, "N": -3, },
        "N": {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, },
    }
    :param extGapPenalty: True if extended gaps should be penalized other way than opened gaps
    :return: scoring matrix
    """
    matrix = [ [0 for x in range(len(s1) + 1) ] for y in range(len(s2) + 1)]

    if extGapPenalty:
        matrix[0][0] = [(0, "l", "")]

        for x in range(1, len(s1) + 1):
            matrix[0][x] = [(GE * x, "l", "g")]
        for y in range(1, len(s2) + 1):
            matrix[y][0] = [(GE * y, "u", "g")]

        for y in range(1, len(s2) + 1):
            for x in range(1, len(s1) + 1):
                matrix[y][x] = allMax((substMatrix[s1[x-1]][s2[y-1]] if substMatrix else (matrix[y-1][x-1][0][0] + M if s2[y-1] == s1[x-1] else matrix[y-1][x-1][0][0] + S), "d", ""),
                                      (matrix[y][x-1][0][0] + G, "l", "g")
                                        if [z[2] if z[1] == "l" else "" for z in matrix[y][x-1]][0] != "g"
                                        else (matrix[y][x-1][0][0] + GE, "l", "g"),
                                      (matrix[y-1][x][0][0] + G, "u", "g")
                                        if [z[2] if z[1] == "u" else "" for z in matrix[y-1][x]][0] != "g"
                                        else (matrix[y-1][x][0][0] + GE, "u", "g"),
                                      key=lambda x: x[0])
    else:
        matrix[0] = [ [(0, "l")] for x in range(len(s1) + 1) ]
        for i in range(len(matrix)):
            matrix[i][0] = [(0, "u")]

        for y in range(1, len(s2) + 1):
            for x in range(1, len(s1) + 1):
                matrix[y][x] = allMax((substMatrix[s1[x-1]][s2[y-1]] if substMatrix else (matrix[y-1][x-1][0][0] + M if s2[y-1] == s1[x-1] else matrix[y-1][x-1][0][0] + S), "d"),
                                      (matrix[y][x-1][0][0] + G, "l"),
                                      (matrix[y-1][x][0][0] + G, "u"),
                                      key=lambda x: x[0])

    return matrix

def needleman(s1, s2, substMatrix=None, extGapPenalty=False, scoringMatrix=None, randomizedBacktrack=False, returnAlternatives=False, returnMatrix=False):
    """
    Returns dictionary containing the two strings alignment or generator containing all the alignments and alignment score.
    Ways in the score matrix can be randomized (i.e. if somewhere on the way are two or more possible following ways, one is randomly chosen).
    Can also return score matrix and you can use your own matrix if needed.

    :param s1: first string to align
    :param s2: second string to align
    :param substMatrix: substitution matrix to use. Expected format:
    {
        "A": {"A": +1, "T": -3, "C": -3, "G": -3, "N": -3, },
        "T": {"A": -3, "T": +1, "C": -3, "G": -3, "N": -3, },
        "C": {"A": -3, "T": -3, "C": +1, "G": -3, "N": -3, },
        "G": {"A": -3, "T": -3, "C": -3, "G": +1, "N": -3, },
        "N": {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, },
    }
    :param extGapPenalty: True if extended gaps should be penalized other way than opened gaps
    :param matrix: scoring matrix
    :param randomizedBacktrack: if True and alternative aligments exist, return random one
    :param returnAlternatives: if True return generator which yields all alternative alignments
    :param returnMatrix: if True return scoring matrix
    :return: {
            "matrix": score matrix if returnMatrix,
            "score": alignment score (i.e. score in right bottom corner of score matrix),
            "alternatives": generator containing all the alignments if returnAlternatives,
            "s1": first aligned string if not returnAlternatives,
            "s2": second aligned string if not returnAlternatives
    }

    >>> n = needleman("TACGT", "CTAGGAT", substMatrix=DNA)
    >>> n["s1"]
    '-TACG-T'
    >>> n["s2"]
    'CTAGGAT'
    >>> n["score"]
    1

    Tested on paper:

    >>> matrix = [
    ... [ [(0, "u")], [(0, "l")], [(0, "l")], [(0, "l")], [(0, "l")],           [(0, "l")] ],
    ... [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "u")],           [(0, "u")] ],
    ... [ [(0, "u")], [(0, "d")], [(0, "d")], [(0, "d")], [(0, "d"), (0, "u")], [(0, "u")] ],
    ... [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "u")],           [(0, "u")] ],
    ... [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "d"), (0, "u")], [(0, "u")] ],
    ... [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "u")],           [(0, "d")] ]]
    >>> s1 = "TACGT"
    >>> s2 = "ACCTG"
    >>> align = needleman("TACGT", "ACCTG", scoringMatrix=matrix)
    >>> align["s1"]
    'TACGT'
    >>> align["s2"]
    'ACCTG'

    >>> align_alts = needleman("TACGT", "ACCTG", scoringMatrix=matrix, returnAlternatives=True)
    >>> for alt in align_alts["alternatives"]:
    ...     alt
    ('TACGT', 'ACCTG')
    ('TACG--T', '--ACCTG')
    ('TACG----T', '----ACCTG')
    """

    if scoringMatrix is None:
        scoringMatrix = generateMatrix(s1, s2, substMatrix=substMatrix, extGapPenalty=extGapPenalty)

    if returnAlternatives:
        return {
            "matrix": scoringMatrix if returnMatrix else None,
            "score": scoringMatrix[-1][-1][0][0],
            "alternatives": _backtrackAlts(scoringMatrix, s1, s2, len(s1), len(s2)),
            "s1": None,
            "s2": None
        }
    else:
        al1, al2 = _backtrackRandomized(scoringMatrix, s1, s2) if randomizedBacktrack else _backtrack(scoringMatrix, s1, s2)
        return {
                "matrix": scoringMatrix if returnMatrix else None,
                "score": scoringMatrix[-1][-1][0][0],
                "alternatives": None,
                "s1": al1,
                "s2": al2
            }

M = 1 # match
S = -1 # substitution (mismatch)
G = -2 # gap open penalty
GE = -1 # gap extend penalty

DNA = {
        "A": {"A": +1, "T": -3, "C": -3, "G": -3, "N": -3, },
        "T": {"A": -3, "T": +1, "C": -3, "G": -3, "N": -3, },
        "C": {"A": -3, "T": -3, "C": +1, "G": -3, "N": -3, },
        "G": {"A": -3, "T": -3, "C": -3, "G": +1, "N": -3, },
        "N": {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, },
}

"""
n = needleman("TACGT", "CTAGGAT", returnMatrix=True, returnAlternatives=True)

print("TACGT", "CTAGGAT")
for x in n["alternatives"]:
    print(x[0])
    print(x[1])
    print()

matrix = [
        [ [(0, "u")], [(0, "l")], [(0, "l")], [(0, "l")], [(0, "l")],           [(0, "l")] ],
        [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "u")],           [(0, "u")] ],
        [ [(0, "u")], [(0, "d")], [(0, "d")], [(0, "d")], [(0, "d"), (0, "u")], [(0, "u")] ],
        [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "u")],           [(0, "u")] ],
        [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "d"), (0, "u")], [(0, "u")] ],
        [ [(0, "u")], [(0, "d")], [(0, "u")], [(0, "d")], [(0, "u")],           [(0, "d")] ]
    ]

print("TACGT", "ACCTG")
align = needleman("TACGT", "ACCTG", matrix=matrix, returnMatrix=True, returnAlternatives=False)
print(align["s1"])
print(align["s2"])
print()
align_alts = needleman("TACGT", "ACCTG", matrix=matrix, returnMatrix=True, returnAlternatives=True)

for i, alt in enumerate(align_alts["alternatives"]):
    print("alt: ", i + 1)
    print(alt[0])
    print(alt[1])
    print()

align2 = needleman("TAAGTAAAAGTCGATTCGGGATCCGTAGCCGTACGTGCAATGCATCGCTAACCAAAGTGCCCGATGCCAATCGT",
                       "ACCTGAGGGTTAATATAGCTCCGCCCTTAAAGATCGCTAGACTCCCGATAGGTCTCAGACT", returnMatrix=True, returnAlternatives=False)

print(align2["s1"])
print(align2["s2"])
print()

align2_alts = needleman("TAAGTAAAAGTCGATTCGGGATCCGTAGCCGTACGTGCAATGCATCGCTAACCAAAGTGCCCGATGCCAATCGT",
                       "ACCTGAGGGTTAATATAGCTCCGCCCTTAAAGATCGCTAGACTCCCGATAGGTCTCAGACT", returnMatrix=True, returnAlternatives=True)

for i, alt in enumerate(align2_alts["alternatives"]):
    print("alt: ", i + 1)
    print(alt[0])
    print(alt[1])
    if i == 0: break
"""

"""import doctest
doctest.testmod()

#generateMatrix("TACGT", "ACCTG", extGapPenalty=True)
print("TACGT", "CTAGGAT")
align = needleman("TACGT", "CTAGGAT", returnMatrix=True)
print(align["score"])
print(align["s1"])
print(align["s2"])
pprint(align["matrix"], width=200)"""

alts = []
s1 = "ATCTAATCCTAAGCCCTACGGTACTCTCAAACTCTGGGCTCAGTTGCTAGCTTGATATAGCTTAGCGCTA"
s2 = "GCATGCGCGAAACTGATACGCTAAGGTCATAGGCATAGCTACAGTCGTCGATACGCTAGCTTGTGATAGTCGCTCGATAGCTGATCGCTAGCGAT"
for x in needleman(s1, s2, returnAlternatives=True)["alternatives"]:
    alts.append(x)
print(alts)
print(len(alts))