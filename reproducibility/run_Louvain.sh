#!/bin/bash

DATASETS=('Zeisel' 'Baron' 'Chen' 'LaManno/Embryo' 'LaManno/ES' 'LaManno/MouseEmbryo')
MY_NORM=('RawCounts' 'TPM' 'DCA' 'Deconvolution' 'MAGIC' 'Sanity' 'SAVER' 'scImpute' 'sctransform' 'scVI' 'SanityErrorbar')

K=(5 10 30 50 100)

for d in ${DATASETS[@]}; do
	echo $d
	for n in ${MY_NORM[@]}; do
		echo -e "\t${n}"
		for k in ${K[@]}; do
			knn_graph=data/$d/${n}_knn${k}.txt
			Louvain_tree=data/$d/${n}_Louvain_tree_K${k}.txt
			Louvain_hierarchy=data/$d/${n}_Louvain_hierarchy_K${k}.txt
			if [ -e $knn_graph ] && [ ! -e $Louvain_tree ]; then
				echo -e -n "\t\t${k} "
				knn_graph_bin='tmp_graph.bin'
				Louvain_convert -i $knn_graph -o $knn_graph_bin
				Louvain_community $knn_graph_bin -l -1 > $Louvain_tree
				Louvain_hierarchy $Louvain_tree > $Louvain_hierarchy
			fi
		done
		echo
	done
done
