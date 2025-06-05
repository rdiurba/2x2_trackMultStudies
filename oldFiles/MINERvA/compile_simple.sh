#compile the CAF plotter


g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` matchDSMinerva.cc -o matchDSMinerva -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` matchUSMinerva.cc -o matchUSMinerva -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict



g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` matchDSPunchMinerva.cc -o matchDSPunchMinerva -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict



g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` matchDSStopMinerva.cc -o matchDSStopMinerva -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict






if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi
