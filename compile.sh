
g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` pandora_sel.cc -o pandora_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` dlp_sel.cc -o dlp_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` minerva_spine.cc -o minerva_spine -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` pandora_minerva.cc -o minerva_pandora -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` minerva_spineSplit.cc -o minerva_spineSplit -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` pandora_minervaSplit.cc -o minerva_pandoraSplit -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi

