
#!/usr/bin/env python

import sys
modules_path = '/home/gencardio/tfm_eli/PIPELINE/modules'
sys.path.append(modules_path)

from  trimming import trimming_fastqs, fastp_metrics
from fastqc import fastqc_trimmed
from mapping import download_reference_genome,sorted_bam_files, mapp_sort_eliminate, preprocessing_variant_calling, eliminate_duplicates, index_bam, download_reference_chromosome
from metrics import  percent_duplications, collect_insert_size_metrics, collect_alignment_summary_metrics, Create_sequence_dictionary,obtain_bed_file, mosdepth, obtain_metrics_csv
from variant_calling import GATK_module, liftovervcf, merge_vcf, merge_vcfs_stats, filter_mutect_calls
from annotation import annotation 


if __name__ == "__main__":


# ---- TRIMMING AND FASTQC modules

    # input_dir_fastqs = "/media/gencardio/KINGSTON/Samples"
    # trimmed_fastqs_dir = "/media/gencardio/KINGSTON/Samples/Results"



    #Trimming fastqs

    trimming_fastq = trimming_fastqs(input_dir_fastqs, trimmed_fastqs_dir)
    
    #Quality report from fastqs

    fastqc = fastqc_trimmed(trimmed_fastqs_dir)

    #fastp_metrics

    fastqs_metrics = fastp_metrics(trimmed_fastqs_dir)

# ---- MAPPING 

    #Download reference genome hg38

    # genome_reference = 'hg38'
    genome_reference = 'hg38'
    download_genome_dir = f"/media/gencardio/KINGSTON/{genome_reference}"
    genome_hg38_path = download_reference_genome(download_genome_dir,genome_reference)
    

    #Obtain the sorted.bam file for each sample and eliminate unecessary files obtained during the mapping step

    final_sorted_bam_files_hg38 = sorted_bam_files(trimmed_fastqs_dir, genome_reference)
    mapping_pipeline = mapp_sort_eliminate(final_sorted_bam_files_hg38, trimming_fastq, genome_hg38_path)


    #Download mithocondrial chromosome
    
    genome_version = 'hg38'
    chromosome_name = 'chrM'
    download_genome_dir = f"/media/gencardio/KINGSTON/{genome_version}/{chromosome_name}"
    chrM_hg38 = download_reference_chromosome(download_genome_dir,genome_version,chromosome_name)

    picard_path = "/home/gencardio/tfm_eli/PIPELINE/picard.jar"

    chrM_dictionary = Create_sequence_dictionary(picard_path, chrM_hg38)

# ---- VARIANT CALLING PRE PROCESING AND METRICS


    perl_script = "/home/gencardio/tfm_eli/PIPELINE/modules/extract_chrM.pl"
    preprocesing_vcf = preprocessing_variant_calling(mapping_pipeline, perl_script)


    #Remove duplicates with picard

    remove_duplicates = eliminate_duplicates(preprocesing_vcf, picard_path)

    #Index sorted bam files
    index_bam_files = index_bam(remove_duplicates)


    #METRICS

    #Duplications %

    percent_duplications_metrics = percent_duplications(index_bam_files, fastqs_metrics)

    #Collect Insert Size metrics

    insert_size_metrics = collect_insert_size_metrics(remove_duplicates, picard_path, percent_duplications_metrics)

    #Collect Alingment Summary Metrics

    alignment_summary_metrics = collect_alignment_summary_metrics(remove_duplicates, picard_path,  genome_hg38_path, insert_size_metrics)

    #Collect mosdepth metrics

    genome_reference_chrM = 'chrM'
    sequence_dictionary = Create_sequence_dictionary(picard_path, chrM_hg38)
    bed_file_path = obtain_bed_file(sequence_dictionary,genome_reference_chrM)
    mosdepth_reports = mosdepth(remove_duplicates, bed_file_path, alignment_summary_metrics)

    #Obtain final metrics .csv

    Metrics_excel = obtain_metrics_csv(mosdepth_reports, trimmed_fastqs_dir)

# ---- VARIANT CALLING

    #GATK 

    gatk_version = '4.6.1.0'

    #Obtain vcf

    reference_path_chrM_shifted = "/media/gencardio/KINGSTON/hg38/hg38_with_shifted_chrM.fa"
    genome_reference_shifted = 'chrM_shifted'
    GATK_module_chrM = GATK_module(remove_duplicates, gatk_version, genome_reference_chrM, genome_reference_shifted, genome_hg38_path, reference_path_chrM_shifted)

    #Obtain the OA flags for the chrM_shifted vcf

    shiftedtohg38 = '/media/gencardio/KINGSTON/hg38/ShiftBack.chain'

    add_OA_flags = liftovervcf(GATK_module_chrM, picard_path, genome_hg38_path, shiftedtohg38)

    #Merge control and non-control region vcf files

    merge_vcf_files = merge_vcf(add_OA_flags, picard_path)

    #Merge control and non-control region vcf stats files

    merge_vcf_stats = merge_vcfs_stats(GATK_module_chrM, gatk_version)

    #Filter Mutect Calls 

    filtered_vcf = filter_mutect_calls(merge_vcf_files, gatk_version, genome_hg38_path)



# ------ ANNOTATION

    #Annotate vcf

    annotation_script = "/home/gencardio/tfm_eli/PIPELINE/modules/Annotation/annotate_mtdna.pl"
    annotated_vcf = annotation(annotation_script, filtered_vcf)

