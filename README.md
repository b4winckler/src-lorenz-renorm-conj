# What is this?

Programs that accompany the article "The Lorenz Renormalization Conjecture".


# Building

Run `make` to build the executables.

Dependencies:

*   A C++ compiler (tested on macOS 10.12.6 with Apple's devtools)
*   The GNU MPFR and GMP libraries (tested with macOS Homebrew installs)

If you are running macOS and do not have these installed, go to the
[Homebrew website](`https://brew.sh`) and follow the installation instructions.
Once Homebrew has finished installing, install mpfr and gmp with
`brew install gmp mpfr`.  Now you should be able to run `make`.


# Running

The program `fixedpt` locates renormalization fixed points; `deriv` calculates
the eigenvalues of the derivative of the renormalization operator; `renorm`
evaluates renormalizations for plotting.


## Example: monotone (2,1)-type

Locate the fixed point of monotone (2,1)-type:

    $ ./fixedpt LRR RL 0.5 2 100 256 100 > fixedpt-2-1.txt

As the program runs it displays errors from the various algorithms; use this to
gauge whether it is converging or not.  Increase the last parameter if it did
not have time to converge, or stop prematurely by typing CTRL-C.

Plot the first two renormalization to a file:

    $ ./renorm LRR RL 1 30 256 < fixedpt-2-1.txt > renorm-2-1.txt

Inspect the eigenvalues:

    $ ./deriv LRR RL 256 < fixedpt-2-1.txt

This should print something like:

    23.1366 12.1264 0.10663

The first two are the unstable eigenvalues associated with moving the critical
values; the last is the eigenvalue associated with moving the critical point.
This is an example of type A in the renormalization conjecture.

Plot column 1 vs 2 (for f), and column 3 vs 4 (for Rf).  For example, here's
how to do it in R:

    > df = read.table('renorm-2-1.txt', header=TRUE)
    > plot(df$x0, df$y0, type='l', col='blue')

You should see a blue graph pop up.  Now add a plot the renormalization:

    > lines(df$x1, df$y1, col='red')

There should now be a red graph obscuring the blue one (since it is a fixed
point the graphs of f and Rf are identical).


## Example: a period-2 orbit of monotone (8,2)-type

Lets try finding a period-2 orbit.  Start with finding the (8,2) fixed point:

    $ ./fixedpt LRRRRRRRR RLL 0.5 2 1000 512 100 > fixedpt-8-2.txt

Now look for the period-2 point:

    $ ./fixedpt LRRRRRRRRRLLRLLRLLRLLRLLRLLRLLRLL \
    > RLLLRRRRRRRRLRRRRRRRR 0.7 2 1000 512 100 > period2-8-2.txt

Note that the critical point is not the same for the two by inspecting the
first column of the second row of the two files just created.  It should be
about 0.14 for the fixed point and 0.75 for the period-2 orbit.

To plot the first two renormalizations of the period-2 point we renormalize
using the (8,2) type; not the twice (8,2) type used to generate the period-2
point:

    $ ./renorm LRRRRRRRR RLL 2 30 512 < period2-8-2.txt > renorm-twice-8-2.txt

Now, plot column 1 vs 2 (for f), column 3 vs 4 (for R(f)), and column 5 vs 6
(for R^2(f)).  For example, in R issue these commands:

    > df = read.table('renorm-twice-8-2.txt', header=TRUE)
    > plot(df$x0, df$y0, type='l', col='blue')
    > lines(df$x1, df$y1, col='red')
    > lines(df$x2, df$y2, col='green')

You should see a green and a red graph; the green will have covered over the
blue since R^2(f) = f.


## Example: scripting

See the `scripts/` folder for scripts on how to automate the generation of
eigenvalues for many (a,b) times at once.


# License

Copyright (c) 2018 Bjorn Winckler, see COPYING for license details.
