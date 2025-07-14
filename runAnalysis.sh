source compile.sh

./pandora_sel Part1MiniRun64.txt testPandoraPart1.root 1
./pandora_sel Part2MiniRun64.txt testPandoraPart2.root 1
./dlp_sel Part1MiniRun64.txt testSPINEPart1.root 1
./dlp_sel Part2MiniRun64.txt testSPINEPart2.root 1

./dlp_sel fhcMCPart1.txt testSPINEFHCPart1.root 1
./dlp_sel fhcMCPart2.txt testSPINEFHCPart2.root 1

./pandora_sel fhcMCPart1.txt testPandoraFHCPart1.root 1
./pandora_sel fhcMCPart2.txt testPandoraFHCPart2.root 1

root -l -b -q plotSPINE.C
root -l -b -q plotSPINEFHC.C

root -l -b -q plotPandora.C
root -l -b -q plotPandoraFHC.C