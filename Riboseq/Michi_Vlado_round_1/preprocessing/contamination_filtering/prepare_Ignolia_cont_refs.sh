#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq


Ignolia_DIR="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper"
Ensembl_cont_DIR="/projects/splitorfs/work/Riboseq/data/contamination"


################################################################################
# RUN GET NCBI rRNA                                                            #
################################################################################
NCBI_cont_fasta="${Ignolia_DIR}"/NCBI_rRNAs.fasta

if [ ! -f ${NCBI_cont_fasta} ]; then
    python retrieve_Ignolia_contamination_sequences.py ${NCBI_cont_fasta}
fi


################################################################################
# CONCATEANTE rRNA reference with mtrRNA from Ensembl                          #
################################################################################
cat ${NCBI_cont_fasta} ${Ensembl_cont_DIR}/Mt_rRNA_trans_110.fasta > ${Ignolia_DIR}/rRNA_ref_NCBI_Ens.fasta


################################################################################
# RUN GET ncRNA Ensembl                                                        #
################################################################################
Rscript /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/preprocessing/contamination_filtering/ncRNA_Ensembl_fasta.R

cat /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ensembl_all_nc_biotypes.fasta \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/gencode.v44.lncRNA_transcripts.fa \
 > /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta