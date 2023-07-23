#!/bin/bash

# Github
../../build/release/tools/dump-graph -i git_web_ml/musae_git_edges.csv --distribution uniform -d LT -o github_LT_binary.txt --scale-factor 0.1 --dump-binary

# HepPh
../../build/release/tools/dump-graph -i cit-HepPh.txt --distribution uniform -d LT -o HepPh_LT_binary.txt --scale-factor 0.1 --dump-binary

# DBLP
../../build/release/tools/dump-graph -i com-dblp.ungraph.txt --distribution uniform -d LT -o DBLP_LT_binary.txt --scale-factor 0.1 --dump-binary -u

# Pokec
../../build/release/tools/dump-graph -i soc-pokec-relationships.txt --distribution uniform -d LT -o Pokec_LT_binary.txt --scale-factor 0.1 --dump-binary

# livejournal
../../build/release/tools/dump-graph -i soc-LiveJournal1.txt --distribution uniform -d LT -o livejournal_LT_binary.txt --scale-factor 0.1 --dump-binary

# orkut-small
../../build/release/tools/dump-graph -i com-orkut.ungraph.txt --distribution uniform -d LT -o orkut_small_LT_binary.txt --scale-factor 0.1 --dump-binary -u

# orkut-big
../../build/release/tools/dump-graph -i download.tsv.orkut-groupmemberships.tar --distribution uniform -d LT -o orkut_big_LT_binary.txt --scale-factor 0.1 --dump-binary -u

# wikipedia
../../build/release/tools/dump-graph -i download.tsv.wikipedia_link_en.tar --distribution uniform -d LT -o wikipedia_LT_binary.txt --scale-factor 0.1 --dump-binary -u

# friendster
../../build/release/tools/dump-graph -i com-friendster.ungraph.txt --distribution uniform -d LT -o friendster_LT_binary.txt --scale-factor 0.1 --dump-binary
