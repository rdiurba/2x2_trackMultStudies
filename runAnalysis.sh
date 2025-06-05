source compile.sh

./pandora_sel Part1MiniRun64.txt testPandoraPart1.root 1
./pandora_sel Part2MiniRun64.txt testPandoraPart2.root 1
./dlp_sel Part1MiniRun64.txt testSPINEPart1.root 1
./dlp_sel Part2MiniRun64.txt testSPINEPart2.root 1
root -l -b -q plotSPINE.C
root -l -b -q plotPandora.C