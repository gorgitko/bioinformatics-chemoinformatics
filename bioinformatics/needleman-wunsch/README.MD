# This is my implementation of Needleman-Wunsch alignment algorithm.

Supports:

* substitution matrix (already contains simple matrix for DNA)
* alternative backtrack ways
* randomized backtrack way when alternative ways exist
* simple constant gap extend penalty (e.g. gap open is for -2 and each gap extend for -1)