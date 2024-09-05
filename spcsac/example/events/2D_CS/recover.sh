#!/bin/bash
NUM=1200
inputdir='../event5/Original_npy/'
outdir='../event5/Recovered_npy/'
filesuffix='npy'
command='python 2DCS_pos.py'
Rsuffix='R.'$filesuffix
Tsuffix='T.'$filesuffix
Zsuffix='Z.'$filesuffix
RTsuffix='RT.'$filesuffix
RZsuffix='RZ.'$filesuffix
RTZsuffix='RTZ.'$filesuffix
pos='../event5/pos.npy'
for i in {200..1200}
do
    inputR=$inputdir$i$Rsuffix
    inputT=$inputdir$i$Tsuffix
    inputZ=$inputdir$i$Zsuffix
    inputRT=$inputdir$i$RTsuffix
    inputRZ=$inputdir$i$RZsuffix
    inputRTZ=$inputdir$i$RTZsuffix
    recover_suffix='Recovered_'
    Rrecovered_outfile=$outdir$i$recover_suffix$Rsuffix
    Trecovered_outfile=$outdir$i$recover_suffix$Tsuffix
    Zrecovered_outfile=$outdir$i$recover_suffix$Zsuffix
    RTrecovered_outfile=$outdir$i$recover_suffix$RTsuffix
    RZrecovered_outfile=$outdir$i$recover_suffix$RZsuffix
    RTZrecovered_outfile=$outdir$i$recover_suffix$RTZsuffix

    $command $inputR $Rrecovered_outfile $pos
    wait
    $command $inputT $Trecovered_outfile $pos
    wait
    $command $inputZ $Zrecovered_outfile $pos
    wait
    $command $inputRT $RTrecovered_outfile $pos
    wait
    $command $inputRZ $RZrecovered_outfile $pos
    wait
    $command $inputRTZ $RTZrecovered_outfile $pos
    wait
    printf "\rProcessing: %4d/%-4d" $i $NUM
done
