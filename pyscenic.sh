#!/bin/sh
#PBS -q tmp
#PBS -l mem=60gb,walltime=500:00:00,nodes=1:ppn=24
#HSCHED -s pop+pyscenic+homo

cd /xtdisk/POP/pyscenic

###EP
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/EP.adjacencies.tsv /scenicdata/EP.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/EP.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/EP.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/EP_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/EP.matrix.tsv /scenicdata/EP_regulons.csv -o /scenicdata/EP_auc_mtx.csv --num_workers 12

###FIB
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/FIB.adjacencies.tsv /scenicdata/FIB.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/FIB.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/FIB.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/FIB_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/FIB.matrix.tsv /scenicdata/FIB_regulons.csv -o /scenicdata/FIB_auc_mtx.csv --num_workers 12


###SMC
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/SMC.adjacencies.tsv /scenicdata/SMC.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/SMC.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/SMC.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/SMC_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/SMC.matrix.tsv /scenicdata/SMC_regulons.csv -o /scenicdata/SMC_auc_mtx.csv --num_workers 12


###MEP
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/MEP.adjacencies.tsv /scenicdata/MEP.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/MEP.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/MEP.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/MEP_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/MEP.matrix.tsv /scenicdata/MEP_regulons.csv -o /scenicdata/MEP_auc_mtx.csv --num_workers 12


###EC
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/EC.adjacencies.tsv /scenicdata/EC.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/EC.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/EC.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/EC_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/EC.matrix.tsv /scenicdata/EC_regulons.csv -o /scenicdata/EC_auc_mtx.csv --num_workers 12


###LEC
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/LEC.adjacencies.tsv /scenicdata/LEC.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/LEC.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/LEC.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/LEC_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/LEC.matrix.tsv /scenicdata/LEC_regulons.csv -o /scenicdata/LEC_auc_mtx.csv --num_workers 12


###MACROPHAGE
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/MACROPHAGE.adjacencies.tsv /scenicdata/MACROPHAGE.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/MACROPHAGE.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/MACROPHAGE.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/MACROPHAGE_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/MACROPHAGE.matrix.tsv /scenicdata/MACROPHAGE_regulons.csv -o /scenicdata/MACROPHAGE_auc_mtx.csv --num_workers 12


###TCELL
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/TCELL.adjacencies.tsv /scenicdata/TCELL.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/TCELL.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/TCELL.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/TCELL_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/TCELL.matrix.tsv /scenicdata/TCELL_regulons.csv -o /scenicdata/TCELL_auc_mtx.csv --num_workers 12


###BCELL
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/BCELL.adjacencies.tsv /scenicdata/BCELL.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/BCELL.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/BCELL.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/BCELL_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/BCELL.matrix.tsv /scenicdata/BCELL_regulons.csv -o /scenicdata/BCELL_auc_mtx.csv --num_workers 12



###PB
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/PB.adjacencies.tsv /scenicdata/PB.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/PB.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/PB.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/PB_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/PB.matrix.tsv /scenicdata/PB_regulons.csv -o /scenicdata/PB_auc_mtx.csv --num_workers 12



###MAST
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic grn -m grnboost2 --num_workers 20 -o /scenicdata/MAST.adjacencies.tsv /scenicdata/MAST.matrix.tsv /scenicdata/allTFs_hg38.txt 
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic ctx /scenicdata/MAST.adjacencies.tsv /scenicdata/hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /scenicdata/MAST.matrix.tsv --mode "dask_multiprocessing" --output /scenicdata/MAST_regulons.csv --num_workers 20
docker run -it --rm -v /xtdisk/POP/pyscenic:/scenicdata docker.io/aertslab/pyscenic:0.9.15 pyscenic aucell /scenicdata/MAST.matrix.tsv /scenicdata/MAST_regulons.csv -o /scenicdata/MAST_auc_mtx.csv --num_workers 12
















