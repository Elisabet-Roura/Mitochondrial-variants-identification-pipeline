#!/usr/bin/env python

import os
import glob
import sys
import subprocess
import logging


#Function to obtain the FastQC quality report 


def fastqc_trimmed(input_dir_trimmed):

    # Obtain all fastq.gz files and it's path
    sample_directories = os.listdir(input_dir_trimmed)

    
    trimmed_fastq_files_list = [] 

    #Fill trimmed_fastq_files_list with all the trimmed fastq paths

    for sample_directory in sample_directories:
        trimmed_fastq_files = glob.glob(f"{input_dir_trimmed}/{sample_directory}/trimming/*.trimmed.fastq.gz")
        trimmed_fastq_files_list.extend(trimmed_fastq_files)
        

    for trimmed_fastq_file in trimmed_fastq_files_list:

        #Obtain name of the trimmed file

        trimmed_file_name = os.path.basename(trimmed_fastq_file)
        sample_name = trimmed_file_name.split('_')[0]
 
    
        trimmed_output_html = trimmed_file_name.replace('.fastq.gz', '_fastqc.html' )
        trimmed_output_zip = trimmed_file_name.replace('.fastq.gz', '_fastqc.zip' )

        
        #Obtain the path of the trimmed file

        trimmed_path = trimmed_fastq_file.split(f'{trimmed_file_name}')[0]
        trimmed_output_html_path = trimmed_path+'/'+ trimmed_output_html
        trimmed_output_zip_path = trimmed_path+'/' + trimmed_output_zip


        #Execute fastqc
        if  os.path.exists(trimmed_fastq_file):

            if not os.path.exists(trimmed_output_html_path):
                if not os.path.exists(trimmed_output_zip_path):
                    cmd = f'fastqc -o {trimmed_path} {trimmed_fastq_file}'
                    print(cmd)
                    p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)     
                    

                    output = p1.stdout.decode("UTF-8")
                    error = p1.stderr.decode("UTF-8")

                    if p1.returncode != 0:
                        msg = f"ERROR: Could not run FastQC for sample {sample_name}"
                        logging.error(msg)
                        logging.error(error)
                        sys.exit(1)

            else:
                if os.path.exists(trimmed_output_zip_path):
                    print(f'The FastQC quality report for {trimmed_file_name} already exists.')

    return trimmed_fastq_files_list






if __name__ == "__main__":
    fastqc_trimmed()
