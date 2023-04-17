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

wget https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz
gunzip soc-pokec-relationships.txt.gz
../../build/release/tools/dump-graph -i soc-pokec-relationships.txt --distribution uniform -d IC -o Pokec_binary.txt --scale-factor 0.1 --dump-binary
rm soc-pokec-relationships.txt

wget https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gunzip soc-LiveJournal1.txt.gz
../../build/release/tools/dump-graph -i soc-LiveJournal1.txt --distribution uniform -d IC -o livejournal_binary.txt --scale-factor 0.1 --dump-binary
rm soc-LiveJournal1.txt

wget https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz
gunzip com-orkut.ungraph.txt.gz
../../build/release/tools/dump-graph -i com-orkut.ungraph.txt --distribution uniform -d IC -o orkut_small_binary.txt --scale-factor 0.1 --dump-binary -u
rm com-orkut.ungraph.txt

wget http://konect.cc/files/download.tsv.orkut-groupmemberships.tar.bz2
bzip2 -d download.tsv.orkut-groupmemberships.tar.bz2
../../build/release/tools/dump-graph -i download.tsv.orkut-groupmemberships.tar --distribution uniform -d IC -o orkut_big_binary.txt --scale-factor 0.1 --dump-binary -u

wget http://konect.cc/files/download.tsv.wikipedia_link_en.tar.bz2
bzip2 -d http://konect.cc/files/download.tsv.wikipedia_link_en.tar.bz2
../../build/release/tools/dump-graph -i download.tsv.wikipedia_link_en.tar --distribution uniform -d IC -o wikipedia_binary.txt --scale-factor 0.1 --dump-binary -u

wget https://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz
gunzip com-friendster.ungraph.txt.gz
../../build/release/tools/dump-graph -i com-friendster.ungraph.txt --distribution uniform -d IC -o friendster_binary.txt --scale-factor 0.1 --dump-binary
rm com-friendster.ungraph.txt.gz