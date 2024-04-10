#compile the CAF plotter
g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` purity_dlp_multiplicity.cc -o purity_plotter -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` eventsel_nominerva_multiplicity.cc -o eventsel_nominerva -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi
