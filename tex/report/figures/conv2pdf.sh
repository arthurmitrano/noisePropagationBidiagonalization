#!/bin/bash
dir=$1
cd $dir

for k in *
do
   epstopdf $k
done
