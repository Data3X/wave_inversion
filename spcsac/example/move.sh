#!/bin/bash
num=2500
for i in {1..$num}
do
    date1='.20240507.PSV.spc'
    date2='.20240507.SH.spc'
    file1=$i$date1
    file2=$i$date2
    cp $file1 $file2 ../sacfile
    cd ../sacfile
    spcsac &> /dev/null
    printf "\rProcessing: %-4d/%-4d" $i $num
    rm *.spc
    cd ..
done