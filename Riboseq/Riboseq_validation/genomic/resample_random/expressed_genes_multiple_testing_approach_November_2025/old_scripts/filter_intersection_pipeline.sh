#----- ----- #

#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq


# Help message:
usage="
Usage: ./filter_intersection_pipeline.sh [-options] [arguments]

where:
-h			show this help
"

# available options for the programm
while getopts 'b:c:e:g:hi:o:s:t:u:d' option; do
  case "$option" in
    b)
        bam="$OPTARG"
        ;;
    c)
        cds_coordinates="$OPTARG"
        ;;
    d)
        dup=true
        ;;
    e)
        ensembl_gtf="$OPTARG"
        ;;
    g)
        genome_fasta="$OPTARG"
        ;;
    i)
        intersection_input="$OPTARG"
        ;;
    o)
        output_star="$OPTARG"
        ;;
    s)
        sample_name="$OPTARG"
        ;;
    t)
        three_primes="$OPTARG"
        ;;  
    u)
        unique_region_dir="$OPTARG"
        ;;  
    h) 
        echo "$usage"
        exit 1
        ;;
   \?) 
        printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done

if [ ! -d "$output_star"/NMD_genome/${sample_name} ];then
    mkdir "$output_star"/NMD_genome/${sample_name}
fi

if [ ! -e  "$output_star"/NMD_genome/${sample_name}/${sample_name}.RNA_Metrics ]; then
    java -jar ~/tools/picard.jar CollectRnaSeqMetrics \
    I=$bam \
    O="$output_star"/NMD_genome/${sample_name}/${sample_name}.RNA_Metrics \
    REF_FLAT=$ensembl_gtf.refflat \
    STRAND=FIRST_READ_TRANSCRIPTION_STRAND
fi

if [ ! -e  "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_intersect_counts_sorted.bed ]; then
    if [ ! -e  "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv ]; then
        htseq-count -f bam --secondary-alignments ignore \
        -c  "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv \
        --supplementary-alignments ignore $bam $ensembl_gtf
    fi

    if [ ! -e "$output_star"/NMD_genome/${sample_name}/3_primes_genomic_merged_numbered_${sample_name}.bed ]; then
        python filter_bed_file_for_expressed_genes_rnanrom.py \
            $three_primes \
            "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv \
            $ensembl_gtf \
            20
    fi

    if [ ! -e  "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_CDS_subtraction_${sample_name}.bed ]; then
        python filter_bed_file_for_expressed_genes_rnanrom.py \
            $unique_region_dir/Unique_DNA_Regions_genomic_CDS_subtraction.bed \
            "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv \
            $ensembl_gtf \
            20
    fi
        
    cds_coordinates_tpm_filtered="${output_star}/NMD_genome/${sample_name}/"$(basename "${cds_coordinates}" .bed)_${sample_name}.bed
    if [ ! -e ${cds_coordinates_tpm_filtered} ]; then
        python filter_bed_file_for_expressed_genes_rnanrom.py \
            "${cds_coordinates}" \
            "${output_star}/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv" \
            $ensembl_gtf \
            20
    fi

    bash riboseq_coverage_3UTRs_vs_CDS_16_12_25.sh -b "${bam}" -c "${cds_coordinates_tpm_filtered}" \
    -s ${sample_name} -t "$output_star"/NMD_genome/${sample_name}/3_primes_genomic_merged_numbered_${sample_name}.bed

    # implement a filter that only conducts the empiricial intersection pipeline
    # if a certain read depth is found

    # this file should contain all reads mapping to mRNA
    # keep sample if the total count ampping to mRNA is larger than a certain number
    if [[ $dup == true ]]; then
        keep_sample=$(python check_seq_depth.py "${output_star}/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv" 7)
        if [ "$keep_sample" = "True" ]; then
            keep_sample=true
        else
            keep_sample=false
        fi
    else
        keep_sample=true
    fi
    
    if [[ $keep_sample == true  ]]; then
        echo $sample_name
        ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
            $intersection_input \
            "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_CDS_subtraction_${sample_name}.bed \
            "$output_star"/NMD_genome/${sample_name}/3_primes_filtered_for_CDS_distribution_${sample_name}_merged.bed \
            "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD \
            "$output_star"/NMD_genome/${sample_name} \
            ${genome_fasta}

        echo "===================       Sample $sample_name intersected"
    fi
fi

