#!/bin/bash

# Exit on error
set -e
RESULTS_PKL="results-SA.pkl"
for SAMPLE in {0..4999}
do
    /usr/bin/python3 GAP_copy.py ERR13577262_raw_string.gfa output load "$RESULTS_PKL" 1 $SAMPLE
    echo "Processing $SAMPLE"
    # mkdir -p fastas2 fastas3 overlaps
    if [ -f "output_${SAMPLE}.gfa" ]; then

        gfatools asm -u output_${SAMPLE}.gfa > ${SAMPLE}_HADOF.gfa

        python3 gfa2unitigs.py ${SAMPLE}_HADOF.gfa > ${SAMPLE}_merged_unitigs.fa


        minimap2 -t 16 ${SAMPLE}_merged_unitigs.fa ERR13577262.fastq > ${SAMPLE}_overlaps.paf
        racon -t 16 ERR13577262.fastq ${SAMPLE}_overlaps.paf ${SAMPLE}_merged_unitigs.fa > ${SAMPLE}_polished1.fasta
        minimap2 -t 16 ${SAMPLE}_polished1.fasta ERR13577262.fastq > ${SAMPLE}_overlaps.paf
        racon -t 16 ERR13577262.fastq ${SAMPLE}_overlaps.paf ${SAMPLE}_polished1.fasta > ${SAMPLE}_polished2.fasta
        minimap2 -t 16 ${SAMPLE}_polished2.fasta ERR13577262.fastq > ${SAMPLE}_overlaps.paf
        racon -t 16 ERR13577262.fastq ${SAMPLE}_overlaps.paf ${SAMPLE}_polished2.fasta > ${SAMPLE}_polished3.fasta
        minimap2 -t 16 ${SAMPLE}_polished3.fasta ERR13577262.fastq > ${SAMPLE}_overlaps.paf
        racon -t 16 ERR13577262.fastq ${SAMPLE}_overlaps.paf ${SAMPLE}_polished3.fasta > ${SAMPLE}_polished_assembly.fasta

        python3 trim_terminal_overlap_reference_free.py \
            ${SAMPLE}_polished_assembly.fasta \
            sequences/${SAMPLE}_final_ref_free_terminal_trimmed.fa \
            --window 700000 \
            --min-overlap 10000 \
            --min-identity 0.99 \
            --terminal-tolerance 2000

    # rm -r fastas2 fastas3 overlaps
        rm ${SAMPLE}_polished1.fasta ${SAMPLE}_polished2.fasta ${SAMPLE}_polished3.fasta ${SAMPLE}_merged_unitigs.fa ${SAMPLE}_HADOF.gfa ${SAMPLE}_overlaps.paf output_${SAMPLE}.gfa
        quast.py -r assembly.fasta ${SAMPLE}_polished_assembly.fasta -o quast_${SAMPLE}
        rm ${SAMPLE}_polished_assembly.fasta
        cp "quast_${SAMPLE}/report.tsv" "results_columns/${SAMPLE}.tsv"
        rm -rf "quast_${SAMPLE}"
    else
        echo "File output_${SAMPLE}.gfa not found, skipping..."
    fi
done
/usr/bin/python3 merge_reports.py results_columns circ_report_SA2.tsv














