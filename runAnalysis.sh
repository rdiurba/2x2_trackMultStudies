source compile.sh

./pandora_sel Part1MiniRun64.txt testPandoraPart1.root 1
./pandora_sel Part2MiniRun64.txt testPandoraPart2.root 1
./dlp_sel Part1MiniRun64.txt testSPINEPart1.root 1
./dlp_sel Part2MiniRun64.txt testSPINEPart2.root 1
./fakedata_dlp_sel Part2MiniRun64.txt testFakeDataHighEReweight.root 1
./fakedata_dlp_selLowE Part2MiniRun64.txt testFakeDataLowEReweight.root 1

./fakedata_pandora_sel Part2MiniRun64.txt testFakeDataPandoraHighEReweight.root 1
./fakedata_pandora_selLowE Part2MiniRun64.txt testFakeDataPandoraLowEReweight.root 1

root -l -b -q plotSPINE.C
#root -l -b -q plotSPINEFHC.C

root -l -b -q plotPandora.C
#root -l -b -q plotPandoraFHC.C

root -l -b -q plotFakeData.C
root -l -b -q plotFakeDataLowE.C

root -l -b -q plotFakeDataPandora.C
root -l -b -q plotFakeDataPandoraLowE.C