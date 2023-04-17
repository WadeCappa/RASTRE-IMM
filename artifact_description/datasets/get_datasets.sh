#!/bin/bash

wget https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gunzip soc-LiveJournal1.txt.gz
../../build/release/tools/dump-graph -i soc-LiveJournal1.txt --distribution uniform -d IC -o livejournal_binary.txt --scale-factor 0.1 --dump-binary
rm soc-LiveJournal1.txt

wget https://snap.stanford.edu/data/cit-HepPh.txt.gz
gunzip cit-HepPh.txt.gz
../../build/release/tools/dump-graph -i cit-HepPh.txt --distribution uniform -d IC -o HepPh_binary.txt --scale-factor 0.1 --dump-binary
rm cit-HepPh.txt

wget https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz
gunzip com-dblp.ungraph.txt
../../build/release/tools/dump-graph -i com-dblp.ungraph.txt --distribution uniform -d IC -o DBLP_binary.txt --scale-factor 0.1 --dump-binary -u
rm com-dblp.ungraph.txt

