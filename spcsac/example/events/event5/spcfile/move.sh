#!/bin/bash
num=2806
for i in {1..2806}
do
    date1='.200710310330A.PSV.spc'
    date2='.200710310330A.SH.spc'
    file1=$i$date1
    file2=$i$date2
    cp $file1 $file2 ../sacfile/
    cd ../sacfile/
    spcsac &> /dev/null
    printf "\rProcessing: %-4d/%-4d" $i $num
    rm *.spc
    cd ../spcfile/
done