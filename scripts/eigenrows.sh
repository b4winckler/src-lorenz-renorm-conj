#!/bin/sh

# usage: eigenrow.sh C B A1 [A2 [A3 ..]]

# Calculate row of eigenvalues monotone (a,b)-type, for a=A1,A2,.. and b=B.
# Use passed in value for C as starting guess; this needs to be close to zero
# if b is much smaller than a!

# Appends output data if file already exists; this is useful if a run failed
# and needs to be restarted.
# For example, we try to generate from (2,2) up to (40,2):
#
#   $ eigenrow.sh 0.5 2 {2..40}
#
# Oops, lets say it failed at a=35 and higher; open output and delete bad rows,
# then:
#
#   $ eigenrow.sh 0.001 2 {35..40}
#
# Now it will only fill in from were it failed.

C=$1
B=$2

fname="data/eigenrow-b$B.dat"
if [ ! -f $fname ]; then
    printf "a\tb\teval0\teval1\teval2\n" > $fname
fi

for A in ${@:3}; do
    w0=$(printf 'L'; printf 'R%.0s' $(eval "echo {1.."$(($A))"}");)
    w1=$(printf 'R'; printf 'L%.0s' $(eval "echo {1.."$(($B))"}");)
    printf "$A\t$B\t" >> $fname
    ./fixedpt $w0 $w1 $C 2 0 512 100 | ./deriv $w0 $w1 512 >> $fname
done
