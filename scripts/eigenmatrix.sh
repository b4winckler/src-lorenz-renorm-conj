#!/bin/sh

# Calculate matrix of eigenvalues for monotone (a,b)-types

fname="data/eigenmatrix.dat"
printf "a\tb\teval0\teval1\teval2\n" > $fname
for i in {1..9}
do
    for j in $(eval echo "{$i..9}")
    do
        w0=$(printf 'L'; printf 'R%.0s' $(eval "echo {1.."$(($j))"}");)
        w1=$(printf 'R'; printf 'L%.0s' $(eval "echo {1.."$(($i))"}");)
        printf "$j\t$i\t" >> $fname
        ./fixedpt $w0 $w1 0.5 2 0 256 100 | ./deriv $w0 $w1 256 >> $fname
    done
done
