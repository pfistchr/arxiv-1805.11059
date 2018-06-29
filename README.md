# Introduction

This is supplementary material for the paper 'Testing Against Independence and a Rényi Information Measure' ([arXiv:1805.11059](https://arxiv.org/abs/1805.11059)) to verify the numerical results in Example 14.

The following software versions were used (make sure you have the development headers installed when you want to compile code):

* [GCC](https://gcc.gnu.org/) 6.3
* [MPFI](http://perso.ens-lyon.fr/nathalie.revol/software.html) 1.5.1
* [MPFR](https://www.mpfr.org/) 3.1.5
* [Mathematica](https://www.wolfram.com/mathematica/) 11.2.0

## Upper Bounds (57) and (58)

To check (57), a specific $R_{XY}$ is shown to satisfy (62) and (63).
Similarly, (58) is established.

### Verification with MPFI

Perform the following steps in the `upperbound` directory:

```
~/arxiv-1805.11059/upperbound$ g++ -O2 -Wall -Wextra -std=c++11 mpfi.cpp -lmpfr -lmpfi -o mpfi
~/arxiv-1805.11059/upperbound$ ./mpfi
0.02973937988281249999
0.81314766615572540082
0.03039550781249999997
0.81022073260010899686
finish
~/arxiv-1805.11059/upperbound$
```

The verification is successful if the program outputs `finish`.
(The expected output is also in the file `mpf.log`.)

### Verification with MPFR

Perform the following steps in the `upperbound` directory:

```
~/arxiv-1805.11059/upperbound$ g++ -O2 -Wall -Wextra -std=c++11 mpfr.cpp -lmpfr -o mpfr
~/arxiv-1805.11059/upperbound$ ./mpfr
0.02973937988281249999
0.81314766615572540082
0.03039550781249999997
0.81022073260010899686
finish
~/arxiv-1805.11059/upperbound$
```

The verification is successful if the program outputs `finish`.
(The expected output is also in the file `mpf.log`.)

### Verification with Mathematica

Perform the following step in the `upperbound` directory:

```
~/arxiv-1805.11059/upperbound$ math -script mathematica.wl
True
True
[...]
True
finish
~/arxiv-1805.11059/upperbound$
```

The verification is successful if all values are `True`, i.e., if the output is identical to the content of the file `mathematica.log`.

## Lower Bound (59)

To check (59), the method described in the paper is used.
A stack initially contains the set $Q$, and the following operations are performed until the stack is empty, i.e., until all subsets of $Q$ have been verified:

* Pop the top element from stack and call it $Q_i$.
* Tighten the $l$'s and $u$'s from (67) and (68), i.e., make sure that for every inequality in the RHS of (67) and (68) there exists a $Q_X Q_Y$ from $Q_i$ that satisfies the inequality with equality.
This operation does not modify the set $Q_i$.
After the operation, the extreme points (see Remark 17) can easily be identified: $Q_X$ is an extreme point if, and only if, one coordinate is equal to the corresponding $l$ and a different coordinate is equal to the corresponding $u$; the analogous statement holds for $Q_Y$.
* Then, either verify with Lemma 16 that $Q_i$ satisfies the lower bound (the `v...` lines in `input.txt` provide the $\alpha$ and $\beta$'s),
* or partition the set $Q_i$ into two subsets and push these subsets to the stack (the `a`...`f` lines in `input.txt` select one of the six possible splits).

### Preparation

Due to GitHub's file size limitations, you first need to create the input file in the `lowerbound` directory:

```
~/arxiv-1805.11059/lowerbound$ cat inputA.txt inputB.txt inputC.txt >input.txt
~/arxiv-1805.11059/lowerbound$
```

### Verification with MPFI

Perform the following steps in the `lowerbound` directory:

```
~/arxiv-1805.11059/lowerbound$ g++ -O2 -Wall -Wextra -std=c++11 mpfi.cpp -lmpfr -lmpfi -o mpfi
~/arxiv-1805.11059/lowerbound$ ./mpfi
0.81628386207460053596
0.87784282021563724813
[...]
0.82851412436105931436
finish
~/arxiv-1805.11059/lowerbound$
```

The verification takes around 30 minutes and is successful if the program outputs `finish`.
(The expected output is also in the file `mpf.log`.)

### Verification with MPFR

Perform the following steps in the `lowerbound` directory:

```
~/arxiv-1805.11059/lowerbound$ g++ -O2 -Wall -Wextra -std=c++11 mpfr.cpp -lmpfr -o mpfr
~/arxiv-1805.11059/lowerbound$ ./mpfr
0.81628386207460053596
0.87784282021563724813
[...]
0.82851412436105931436
finish
~/arxiv-1805.11059/lowerbound$
```

The verification takes around 15 minutes and is successful if the program outputs `finish`.
(The expected output is also in the file `mpf.log`.)

### Verification with Mathematica

Not possible at the moment because of a bug in Mathematica (CASE:4034970).
For example, the output in the following test case should be `False` and not `True`:

```
~/arxiv-1805.11059/lowerbound$ math
Mathematica 11.2.0 Kernel for Linux x86 (64-bit)
Copyright 1988-2017 Wolfram Research, Inc.

In[1]:= 10*29^(7151557967478901/9007199254740992) > 145

Out[1]= True

In[2]:=
~/arxiv-1805.11059/lowerbound$
```
