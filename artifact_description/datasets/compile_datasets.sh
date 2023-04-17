#!/bin/bash

../../build/release/tools/dump-graph -i soc-LiveJournal1.txt --distribution uniform -d IC -o livejournal_binary.txt --scale-factor 0.1 --dump-binary

../../build/release/tools/dump-graph -i cit-HepPh.txt --distribution uniform -d IC -o HepPh_binary.txt --scale-factor 0.1 --dump-binary

../../build/release/tools/dump-graph -i com-dblp.ungraph.txt --distribution uniform -d IC -o DBLP_binary.txt --scale-factor 0.1 --dump-binary -u

../../build/release/tools/dump-graph -i soc-pokec-relationships.txt --distribution uniform -d IC -o Pokec_binary.txt --scale-factor 0.1 --dump-binary

../../build/release/tools/dump-graph -i soc-LiveJournal1.txt --distribution uniform -d IC -o livejournal_binary.txt --scale-factor 0.1 --dump-binary

../../build/release/tools/dump-graph -i com-orkut.ungraph.txt --distribution uniform -d IC -o orkut_small_binary.txt --scale-factor 0.1 --dump-binary -u

../../build/release/tools/dump-graph -i download.tsv.orkut-groupmemberships.tar --distribution uniform -d IC -o orkut_big_binary.txt --scale-factor 0.1 --dump-binary -u

../../build/release/tools/dump-graph -i download.tsv.wikipedia_link_en.tar --distribution uniform -d IC -o wikipedia_binary.txt --scale-factor 0.1 --dump-binary -u

../../build/release/tools/dump-graph -i com-friendster.ungraph.txt --distribution uniform -d IC -o friendster_binary.txt --scale-factor 0.1 --dump-binary