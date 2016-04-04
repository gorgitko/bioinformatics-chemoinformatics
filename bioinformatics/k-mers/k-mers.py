import sys
import timeit
from itertools import product
sys.setrecursionlimit(100000)

def generate_kmers(chars, length, _string=""):
    if len(_string) == length:
        yield _string
    else:
        for char in chars:
            yield from generate_kmers(chars, length, _string=_string + char)

for i, x in enumerate(generate_kmers(["A", "T", "C", "G"], 5)):
    print(i, x)

# this does the same
"""
for i, x in enumerate(product("ATGC", repeat=5)):
    print(i, "".join(x))
"""

# time comparison of my implementation vs. standard itertools
# itertools are probably written in C :)

"""
# measure the execution time of generate_kmers()
setup = "from __main__ import kmer"

print(min(timeit.Timer("[i for i in kmer('ATGC', 10)]", setup=setup).repeat(5, 10)))

# measure the execution time of standard itertools.product()
setup = "from itertools import product"

print(min(timeit.Timer("[i for i in product('ATGC', repeat=10)]", setup=setup).repeat(5, 10)))
"""
