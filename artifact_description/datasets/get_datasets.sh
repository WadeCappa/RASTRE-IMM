#!/bin/bash

# Github
wget https://snap.stanford.edu/data/git_web_ml.zip
unzip git_web_ml.zip

# HepPh
wget https://snap.stanford.edu/data/cit-HepPh.txt.gz
gunzip cit-HepPh.txt.gz

# DBLP
wget https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz
gunzip com-dblp.ungraph.txt

# Pokec
wget https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz
gunzip soc-pokec-relationships.txt.gz

# livejournal
wget https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gunzip soc-LiveJournal1.txt.gz

# orkut-small
wget https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz
gunzip com-orkut.ungraph.txt.gz

# orkut-big
wget http://konect.cc/files/download.tsv.orkut-groupmemberships.tar.bz2
bzip2 -d download.tsv.orkut-groupmemberships.tar.bz2

# wikipedia
wget http://konect.cc/files/download.tsv.wikipedia_link_en.tar.bz2
bzip2 -d http://konect.cc/files/download.tsv.wikipedia_link_en.tar.bz2

# friendster
wget https://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz
gunzip com-friendster.ungraph.txt.gz