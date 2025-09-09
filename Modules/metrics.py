#!/usr/bin/env python

import sys
import pandas as pd
import os
import subprocess
import logging
import gzip


#All this module is designed forward to obtain diferent metrics

#Is only performed in the files mapped to the hg38 and not the shifted version


#Function to obtain the % of duplications

def percent_duplications(index_bam_files_list, summary_qc_dict):

    for index_bam_file in index_bam_files_list:

        if not 'shifted' in index_bam_file:

            duplication_file = index_bam_file.replace('_merged.rd.bam.bai','_merged.metrics.txt')
            sample_name = os.path.basename(index_bam_file).split('.')[0]

            found = False
            with open(duplication_file, 'rt', encoding='utf-8', errors='replace') as file: 
                for line in file:
                    if "PERCENT_DUPLICATION" in line:
                        found = True
                        continue
                    if found:

                        tmp = line.strip().split('\t')[-2].replace(',','.')

                        if sample_name not in summary_qc_dict:
                            print(f'ERROR: Sample {sample_name} was not found in the metrics dictionary.')
                            sys.exit(1)

                        summary_qc_dict[sample_name]['duplications_%'] = float(tmp)*100
                        break

    print(summary_qc_dict)
                    

    return summary_qc_dict


#Function to obtain the mean insert size and the standard deviation

def collect_insert_size_metrics(removed_duplicate_files_list, picard_path, metrics_dictionary):

    for removed_duplicate_file in removed_duplicate_files_list:

        if not 'shifted' in removed_duplicate_file:

            metrics_file_name = removed_duplicate_file.replace('.rd.bam','.insert_size_metrics.txt')
            pdf_histogram_file_name =  removed_duplicate_file.replace('.rd.bam','.histogram.pdf')
            sample_name = os.path.basename(removed_duplicate_file).split('.')[0]



            if not os.path.exists(metrics_file_name):

                cmd = f'java -jar {picard_path} CollectInsertSizeMetrics \
                -I {removed_duplicate_file}\
                -O {metrics_file_name}\
                -H {pdf_histogram_file_name}' 
                print(cmd)  
                p_insert = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                

                output_insert = p_insert.stdout.decode("UTF-8")
                error_insert = p_insert.stderr.decode("UTF-8")

                if p_insert.returncode != 0:
                    msg = f"ERROR: Could not run CollectInsertSizeMetrics for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_insert)
                    sys.exit(1)

            else:
                if os.path.exists(metrics_file_name):
                    print(f'CollectInsertSizeMetrics report already exists for {removed_duplicate_file}.')

            found = False
            with open(metrics_file_name, 'r', encoding='utf-8', errors='replace') as file: 
                for line in file:
                    if "MEAN_INSERT_SIZE" in line:
                        found = True
                        continue
                    if found:

                        mean_insert_size = line.strip().split('\t')[5].replace(',','.')
                        sd = line.strip().split('\t')[6].replace(',','.')

                        metrics_dictionary[sample_name]['mean_insert_size'] = mean_insert_size
                        metrics_dictionary[sample_name]['sd']= sd

                        if sample_name not in metrics_dictionary:
                            print(f'ERROR: Sample {sample_name} was not found in the metrics dictionary.')
                            sys.exit(1)

                        break

    print(metrics_dictionary)

    return metrics_dictionary


#Function to obtain the the following metrics: total mapped reads, % mapped reads, pf reads, mean read lenght

def collect_alignment_summary_metrics(removed_duplicate_files_list, picard_path, genome_path, metrics_dictionary):


    for removed_duplicate_file in removed_duplicate_files_list:

        if not 'shifted' in removed_duplicate_file:



            metrics_file_name = removed_duplicate_file.replace('.rd.bam','.aligment_summary_metrics.txt')
            sample_name = os.path.basename(removed_duplicate_file).split('.')[0]

            filtered_bam = removed_duplicate_file.replace('.rd.bam', '.chrM_only.bam')
            filtered_bam_index = filtered_bam + ".bai"
            metrics_file_name = removed_duplicate_file.replace('.rd.bam','.aligment_summary_metrics.txt')

            if not os.path.exists(filtered_bam):

                cmd_filter = f"samtools view -b {removed_duplicate_file} chrM > {filtered_bam}"
                print(f"Filtering BAM for chrM: {cmd_filter}")
                subprocess.run(cmd_filter, shell=True, check=True)

                
                cmd_index = f"samtools index {filtered_bam}"
                print(f"Indexing filtered BAM: {cmd_index}")
                subprocess.run(cmd_index, shell=True, check=True)

            if not os.path.exists(metrics_file_name):

                cmd = f'java -jar {picard_path} CollectAlignmentSummaryMetrics\
                -R {genome_path}\
                -I {filtered_bam}\
                -O {metrics_file_name}'
                print(cmd)
                p_summary = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                

                output_summary = p_summary.stdout.decode("UTF-8")
                error_summary = p_summary.stderr.decode("UTF-8")

                if p_summary.returncode != 0:
                    msg = f"ERROR: Could not run CollectAlignmentSummaryMetrics for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_summary)
                    sys.exit(1)
            else:

                if os.path.exists(metrics_file_name):
                    print(f'CollectAlignmentSummaryMetrics report already exists for {removed_duplicate_file}.')

            found = False
            with open(metrics_file_name, 'r',  encoding='utf-8', errors='replace') as file: 
                for line in file:
                    if "SECOND_OF_PAIR" in line:
                        found = True
                        continue
                    if found:

                        total_reads = line.strip().split('\t')[1].replace(',','.')
                        pf_reads = line.strip().split('\t')[2].replace(',','.')
                        mean_read_lenght = line.strip().split('\t')[15].replace(',','.')

                        inicial_reads = metrics_dictionary[sample_name]['inicial_reads_(M)'] 

                        # metrics_dictionary[sample_name]['total_mapped_reads_(M)'] = float(total_reads)/1000000
                        metrics_dictionary[sample_name]['total_mapped_reads'] = total_reads
                        metrics_dictionary[sample_name][' %_mapped'] = (float(total_reads)/1000000)*100/float(inicial_reads)
                        metrics_dictionary[sample_name]['pf_reads'] = pf_reads
                        # metrics_dictionary[sample_name]['pf_reads_(M)'] = float(pf_reads)/1000000
                        metrics_dictionary[sample_name]['mean_read_lenght'] = mean_read_lenght

                        if sample_name not in metrics_dictionary:
                            print(f'ERROR: Sample {sample_name} was not found in the metrics dictionary.')
                            sys.exit(1)

                        break

    print(metrics_dictionary)

    return metrics_dictionary

def Create_sequence_dictionary(picard_path,genome_path):

    sequence_dictionary = genome_path.replace('.fa','.dict')

    if not os.path.exists(sequence_dictionary):
        cmd = f'java -jar {picard_path} CreateSequenceDictionary -R {genome_path} -O {sequence_dictionary}'
        subprocess.run(cmd, shell= True)
        print(cmd)
    
    index_reference_file = genome_path.replace('.fa','.fa.fai')


    if not os.path.exists(index_reference_file):
        cmd = f'samtools faidx {genome_path}'
        print(cmd)

        p_faidx = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        output_faidx = p_faidx.stdout.decode("UTF-8")
        error_faidx = p_faidx.stderr.decode("UTF-8")

        if p_faidx.returncode != 0:
            msg = f"ERROR: Could not index reference genome with samtools faidx for {genome_path}"
            logging.error(msg)
            logging.error(error_faidx)
            sys.exit(1)

    return sequence_dictionary


def obtain_bed_file (sequence_dictionary, genome_reference):
    
    fai_path = sequence_dictionary.replace('dict', 'fa.fai')
    bed_path = sequence_dictionary.replace('dict',  'bed')

    if not os.path.exists(fai_path):
        msg = f"ERROR: FAI file not found: {fai_path}"
        logging.error(msg)
        sys.exit(1)

    with open(fai_path, 'r') as fai, open(bed_path, 'w') as bed:
        for line in fai:
            fields = line.strip().split('\t')
            if fields[0] == genome_reference:
                bed.write(f"{fields[0]}\t0\t{fields[1]}\n")
                print(f"[INFO] BED file for {genome_reference} created at {bed_path}")
                return bed_path

    msg = f"ERROR: Chromosome '{genome_reference}' not found in {fai_path}"
    logging.error(msg)
    sys.exit(1)


#Function to obtain the the following metrics: mean target coverage, PCT Target Bses x1,2,10,20,30,40,50,100   

def mosdepth(removed_duplicates_files_list, bed_file_path, metrics_dictionary):

    thresholds = [1, 2, 10, 20, 30, 40, 50, 100]

    for removed_duplicate_file in removed_duplicates_files_list:

        if not 'shifted' in removed_duplicate_file:


            sample_name = os.path.basename(removed_duplicate_file).split('.')[0]
            file_path = os.path.dirname(removed_duplicate_file)
            output_prefix = os.path.join(file_path, sample_name)

            # mosdepth_report = output_prefix + '.mosdepth.global.dist.txt'
            region_report = output_prefix + '.regions.bed.gz'

            per_base_report = output_prefix + '.per-base.bed.gz'
            # print(region_report, per_base_report)


            if not os.path.exists(region_report) or not os.path.exists(per_base_report):

                cmd = f'mosdepth --by {bed_file_path} {output_prefix} {removed_duplicate_file}'
                print(cmd)
                p_mosdepth = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    

                output_mosdepth = p_mosdepth.stdout.decode("UTF-8")
                error_mosdepth = p_mosdepth.stderr.decode("UTF-8")

                if p_mosdepth.returncode != 0:
                    msg = f"ERROR: Could not run Mosdepth for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_mosdepth)
                    

            else:
                print(f'Mosdepth report already exists for {sample_name}.')
            
            if os.path.exists(region_report):


                with gzip.open(region_report, 'rt',  encoding='utf-8', errors='replace') as file: 

                    for line in file:

                        if line.split('\t', 1)[0] == 'chrM':
                        
                            mean_chrM_coverage = line.strip().split('\t')[3]

                            metrics_dictionary[sample_name]['mean_target_coverage'] = mean_chrM_coverage

                            
                        
                        else: 
                            print(f'Could not find chrM mean coverage for {sample_name}')
                

            
            if os.path.exists(per_base_report):

                total_bases = 0
                total_coverage = 0
                covered_bases = {t: 0 for t in thresholds}

                with gzip.open(per_base_report, 'rt') as file:
                    for line in file:
                        if not line.startswith("chrM"):
                            continue
                        chrom, start, end, depth = line.strip().split('\t')
                        span = int(end) - int(start)
                        cov = float(depth)

                        total_bases += span
                        total_coverage += cov * span

                        for t in thresholds:
                            if cov >= t:
                                covered_bases[t] += span

                if total_bases == 0:
                    logging.warning(f"[WARN] No hi ha bases per calcular call rate per {sample_name}")
                    continue

                for t in thresholds:
                    rate = 100 * covered_bases[t] / total_bases
                    metrics_dictionary[sample_name][f'Call_rate x{t}'] = rate

                    # print(rate)
                    
    print(metrics_dictionary)
    return metrics_dictionary


# Function to convert the dictionary with all the metrics collected in the diferents functions into a .csv file

def obtain_metrics_csv(metrics_dictionary, csv_path):

    metrics_file = csv_path  + '/Metrics.csv'


    df = pd.DataFrame.from_dict(metrics_dictionary, orient='index')
    print(df)
    df.to_csv(metrics_file, index=True, decimal='.')

    print(metrics_file)
    


if __name__ == "__main__":

    percent_duplications()
    collect_insert_size_metrics()
    collect_alignment_summary_metrics()
    Create_sequence_dictionary()
    obtain_bed_file()
    mosdepth()
    obtain_metrics_csv()
    
