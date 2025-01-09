#!/bin/bash

for i in {2..88..2}
do
 matmult_c.gcc blk $1 $2 $3 $i
done
