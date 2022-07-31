# a179009
Program to count maximally refined partition of a number N (see https://oeis.org/A179009)

The program supports multicore on a single node, but does not support a cluster system on multiple nodes.


```Usage: a179009 [N | -e N | N M]

 - no arguments: batch mode, check sequences in stdin
 - one arguments: computes a(N)
 - two arguments (N<=M): computes a(N)...a(M)
 - '-e N' enumerates non refinable partitions of N


See: 
[https://oeis.org/A179009]
[https://oeis.org/A179009/b179009.txt]
