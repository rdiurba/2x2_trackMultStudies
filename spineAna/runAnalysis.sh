source compile.sh
rm output.root
#./pandora_sel Part1MiniRun65.txt testPandoraPart1.root 1
#./pandora_sel Part2MiniRun65.txt testPandoraPart2.root 1
./dlp_sel Part1MiniRun65.txt testSPINEPart1.root 1
./dlp_fluxRw_sel Part1MiniRun65.txt testFluxRWSPINE.root 1
./dlp_genieRw_sel Part1MiniRun65.txt testGENIERWSPINE.root 1

./dlp_sel Part2MiniRun65.txt testSPINEPart2.root 1
#./fakedata_dlp_sel Part2MiniRun65.txt testFakeDataHighEReweight.root 1
#./fakedata_dlp_selLowE Part2MiniRun65.txt testFakeDataLowEReweight.root 1

root -l -b -q unfoldGENIERW.C
root -l -b -q unfoldFluxRW.C 
root -l -b -q testEff.C
root -l -b -q plotSPINE.C


