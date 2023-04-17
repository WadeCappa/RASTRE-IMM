#!/bin/bash

wget https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gunzip soc-LiveJournal1.txt.gz

wget https://snap.stanford.edu/data/cit-HepPh.txt.gz
gunzip cit-HepPh.txt.gz

wget https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz
gunzip com-dblp.ungraph.txt

wget https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz
gunzip soc-pokec-relationships.txt.gz

wget https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gunzip soc-LiveJournal1.txt.gz

wget https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz
gunzip com-orkut.ungraph.txt.gz

wget http://konect.cc/files/download.tsv.orkut-groupmemberships.tar.bz2
bzip2 -d download.tsv.orkut-groupmemberships.tar.bz2

wget http://konect.cc/files/download.tsv.wikipedia_link_en.tar.bz2
bzip2 -d http://konect.cc/files/download.tsv.wikipedia_link_en.tar.bz2

wget https://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz
gunzip com-friendster.ungraph.txt.gz