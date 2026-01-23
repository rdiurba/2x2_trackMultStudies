source compile.sh

./minerva_pandoraSplit MiniRun64.txt testMinervaPandoraSplit.root 1
./minerva_pandora MiniRun64.txt testMinervaPandora.root 1

./minerva_spineSplit MiniRun64.txt testMinervaSPINESplit.root 1
./minerva_spine MiniRun64.txt testMinervaSPINE.root 1

./minerva_pandoraSplit Sandbox.txt testMinervaPandoraDataSplit.root 0
./minerva_pandora Sandbox.txt testMinervaPandoraData.root 0

./minerva_spineSplit Sandbox.txt testMinervaSPINEDataSplit.root 0
./minerva_spine Sandbox.txt testMinervaSPINEData.root 0

root -l -b -q plotSPINEMinervaSplit.C
root -l -b -q plotSPINEMinerva.C
root -l -b -q plotMinervaData.C
root -l -b -q plotMinervaDataSplit.C
