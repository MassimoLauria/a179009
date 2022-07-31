# a179009
Program to count maximally refined partition of a number N (see https://oeis.org/A179009)

The program supports multicore on a single node, but does not support a cluster system on multiple nodes.

## Compilation
The program is a single C++ file which uses OpenMP. Compile it with something like

```sh
g++ -O3 a179009.cc -fopenmp -o a179009
```

## Usage

By running `a179009 -h` you get follwing help.

```
Usage: a179009 [N | -e N | N M]

 - no arguments: batch mode, check sequences in stdin
 - one arguments: computes a(N)
 - two arguments (N<=M): computes a(N)...a(M)
 - '-e N' enumerates non refinable partitions of N
```

See: 
- https://oeis.org/A179009
- https://oeis.org/A179009/b179009.txt
