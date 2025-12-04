#!/bin/bash
# ChIP-seq analysis pipeline - Peak calling module
# Original author: gxu@lawlab
# v4 changes: removed the quotation marks from the peak output files; removed the summary step
# Adapted for container environment by Salk HPI
set -euo pipefail

# Usage:
#   ./ChIProfiler_M03_CallPeaks.sh <IP_name> <control_name> <tag_suffix> <identifier>
# Arguments:
#   <IP>            : Name of IP tag directory prefix (e.g., CLSY3)
#   <Control>       : Name of control tag directory prefix (e.g., WT)
#   <Suffix>        : Shared suffix for tag directories (e.g., _TagDir)
#   <identifier>    : Identifier used in output filenames (e.g., z07)

################################# parameters
#################################
IP="$1$3"
control="$2$3"
identifier=$4

################################# call peaks with different thresholds
#################################
outputlinker="_vs"
outputsuffix=".txt"

findPeaks $IP -i $control -style factor -region -L 2 -F 2 -center -o $identifier\_testPeaks1_$1$outputlinker$2$outputsuffix &
findPeaks $IP -i $control -style factor -region -L 2.5 -F 2.5 -center -o $identifier\_testPeaks2_$1$outputlinker$2$outputsuffix &
findPeaks $IP -i $control -style factor -region -L 3 -F 3 -center -o $identifier\_testPeaks3_$1$outputlinker$2$outputsuffix &
findPeaks $IP -i $control -style factor -region -L 3.5 -F 3.5 -center -o $identifier\_testPeaks4_$1$outputlinker$2$outputsuffix &
findPeaks $IP -i $control -style factor -region -L 4 -F 4 -center -o $identifier\_testPeaks5_$1$outputlinker$2$outputsuffix
wait

for f in $identifier\_testPeaks*; do grep -v "#" $f | awk '{print $2, $3,$4,$1}' | sed 's/ /\t/g'  > ${f/.txt/.bed}; done
