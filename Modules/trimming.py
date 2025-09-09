#!/usr/bin/env python

import sys
import os
import subprocess
import glob
from os.path import join
import json
import csv
import logging


#Define the function for trimming fastqs

def trimming_fastqs(input_dir, output_dir):

    paired_dictionary = {}

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)


    # Obtain all fastq.gz files and it's path

    fastq_files = glob.glob(f"{input_dir}/*.fastq.gz")

    #Obtain the dictionary of the paired runs for each sample

    sample_list = []
   
    for fastq in fastq_files:
        fastq_name = os.path.basename(fastq)
        sample_name =  fastq_name.split('_')[0]
        if sample_name not in sample_list:
            sample_list.append(sample_name)

        if not sample_name in paired_dictionary:
            paired_dictionary[sample_name] = {}
        if "R1" in fastq:
            paired_dictionary[sample_name]['R1'] = fastq
        if "R2" in fastq:
            paired_dictionary[sample_name]['R2'] = fastq



    for sample_name in paired_dictionary:

        #Create a directory for each sample

        sample_dir = os.path.join(output_dir, sample_name)
        if not os.path.exists(sample_dir):
            sample_dir_trimming = sample_dir +'/trimming'
            os.makedirs(sample_dir_trimming)

        if os.path.exists(sample_dir) and not os.path.exists(f'{sample_dir}/trimming'):
            sample_dir_trimming = sample_dir +'/trimming'
            os.makedirs(sample_dir_trimming)

            

        #fill paired_dictionary 

        if not "R1" in paired_dictionary[sample_name]:
                continue
        if not "R2" in paired_dictionary[sample_name]:
                continue
        

        # Obtain fq1 and fq2 complete path

        fq1 = paired_dictionary[sample_name]['R1'] 
        fq2 = paired_dictionary[sample_name]['R2'] 
        trimmed_fq1_name = os.path.basename(fq1.replace(".fastq.gz", ".trimmed.fastq.gz"))
        trimmed_fq2_name = os.path.basename(fq2.replace(".fastq.gz", ".trimmed.fastq.gz"))


        #Obtain the name for the trimmed files

        paired_dictionary[sample_name]['TRIMMED_R1'] = sample_dir + '/trimming/' + trimmed_fq1_name
        paired_dictionary[sample_name]['TRIMMED_R2'] = sample_dir + '/trimming/' + trimmed_fq2_name


        trimmed_fq1_path = sample_dir +'/trimming/'+trimmed_fq1_name 
        trimmed_fq2_path = sample_dir+'/trimming/'+trimmed_fq2_name
        reports_path = sample_dir+'/trimming/'+sample_name
        

        #cmd to execute fastp

        if not os.path.exists(trimmed_fq1_path) and not os.path.exists(trimmed_fq2_path):
            cmd = f'fastp -i {fq1}\
            -I {fq2}\
            -o {trimmed_fq1_path}\
            -O {trimmed_fq2_path}\
            -h {reports_path}.fastp.html\
            -j {reports_path}.fastp.json'
            print(cmd)
            p_fastp = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output = p_fastp.stdout.decode("UTF-8")
            error = p_fastp.stderr.decode("UTF-8")

            if p_fastp.returncode != 0:
                msg = f" ERROR: Could not run fastp for sample {sample_name}"
                logging.error(msg)
                logging.error(error)
                sys.exit(1)
         
               
    print(paired_dictionary)
    return paired_dictionary



#Obtain dictionary with fastp's metrics


def fastp_metrics(input_dir_json):

    summary_qc_dict = {}
    json_directories = os.listdir(input_dir_json)


    json_files_list = [] 

    #Fill jason_files_list 

    for json_directory in json_directories:
        json_files = glob.glob(f"{input_dir_json}/{json_directory}/trimming/*.fastp.json")
        json_files_list.extend(json_files)
            


    for json_files in json_files_list :

        #Obtain name of the json file and it's path

        json_file =  json_files.split('/')[-1]
        sample_name = json_file.split('.')[0]


        json_path = json_files.split(f'{json_file}')[0]

        os.chdir(json_path)

        with open (json_file) as archivo_json:
            archivo_dict = json.load(archivo_json)

            
            #metrics: q20, q30, total_trimmed_bases, %trimmed_bades, reads_with_adapter, %reads_with_adapter, gc_content

            q20 = archivo_dict['summary']['after_filtering']['q20_rate'] * 100
            q30 = archivo_dict['summary']['after_filtering']['q30_rate'] * 100
            bases_before_filtering = archivo_dict['summary']['before_filtering']['total_bases']
            bases_after_filtering = archivo_dict['summary']['after_filtering']['total_bases']
            total_trimmed_bases = bases_before_filtering - bases_after_filtering 
            trimmed_bases_perc = total_trimmed_bases*100/bases_before_filtering
            total_inicial_reads = archivo_dict['summary']['before_filtering']['total_reads']
            reads_with_adapter = archivo_dict['adapter_cutting']['adapter_trimmed_reads']
            reads_with_adapter_perc = reads_with_adapter*100/ archivo_dict['summary']['before_filtering']['total_reads']
            gc_content = archivo_dict['summary']['after_filtering']['gc_content']

            
            if not sample_name in summary_qc_dict:
                summary_qc_dict[sample_name] = {}

            if not "q20_rate" in summary_qc_dict:
                summary_qc_dict[sample_name]['q20_rate'] = q20
            if not "q30_rate" in summary_qc_dict:
                summary_qc_dict[sample_name]['q30_rate'] = q30
            if not "total_bases_(M)" in summary_qc_dict:
                summary_qc_dict[sample_name]['total_bases_(M)'] = float(bases_before_filtering)/1000000
            if not "total_trimmed_bases_(M)" in summary_qc_dict:
                summary_qc_dict[sample_name]['total_trimmed_bases_(M)'] = float(total_trimmed_bases)/1000000
            if not "trimmed_bases_perc" in summary_qc_dict:
                summary_qc_dict[sample_name]['trimmed_bases_%'] = trimmed_bases_perc
            if not "inicial_reads_(M)" in summary_qc_dict:
                summary_qc_dict[sample_name]['inicial_reads_(M)'] = float(total_inicial_reads)/1000000
            if not "reads_with_adapter" in summary_qc_dict:
                summary_qc_dict[sample_name]['reads_with_adapter'] = reads_with_adapter
            if not "reads_with_adapter_%" in summary_qc_dict:
                summary_qc_dict[sample_name]['read_with_adapter_%'] = reads_with_adapter_perc
            if not "gc_content_%" in summary_qc_dict:
                summary_qc_dict[sample_name]['gc_content_%'] = float(gc_content)*100
                        


    os.chdir(input_dir_json)



    with open('summary_qc_dict.csv', "w", newline="") as f:
        w = csv.DictWriter(f, summary_qc_dict.keys())
        w.writeheader()
        w.writerow(summary_qc_dict)

    f.close()


    return summary_qc_dict


#As a module to import:

if __name__ == "__main__":
    trimming_fastqs()
    fastp_metrics()
    


