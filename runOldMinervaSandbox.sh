source compile.sh

./minerva_pandoraSplit MiniRun64.txt testMinervaPandoraSplit.root 1
./minerva_pandora MiniRun64.txt testMinervaPandora.root 1

./minerva_spineSplit MiniRun64.txt testMinervaSPINESplit.root 1
./minerva_spine MiniRun64.txt testMinervaSPINE.root 1

./minerva_pandoraSplit Sandboxv6.txt testMinervaPandoraDataSplitOld.root 0
./minerva_pandora Sandboxv6.txt testMinervaPandoraDataOld.root 0

./minerva_spineSplit Sandboxv6.txt testMinervaSPINEDataSplitOld.root 0
./minerva_spine Sandboxv6.txt testMinervaSPINEDataOld.root 0

root -l -b -q plotSPINEMinervaSplit.C
root -l -b -q plotSPINEMinerva.C
root -l -b -q plotMinervaData.C
root -l -b -q plotMinervaDataSplit.C
