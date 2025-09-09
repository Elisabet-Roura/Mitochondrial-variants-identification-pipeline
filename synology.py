#!/usr/bin/env python

import sys
import os
import ftputil
from os.path import join
import json
import pandas as pd
import re



def fastqs_identification(path, fastq_dict):

    with ftputil.FTPHost("x", "x", "x") as ftp_host: # access to synology
    # all the runs
        runs = ftp_host.listdir(path)
        # print(runs)
        # obtain the path for all the runs
        for run in runs:
            run_path = path + "/" + run

            if ftp_host.path.isdir(run_path): 
                fastqs = ftp_host.listdir(run_path)
        
                for fastq in fastqs:
                    fastq_name =  fastq.split('_')[0]
                    fastq_path = join(run_path,fastq)  
            
                    if not fastq_name in fastq_dict:
                        fastq_dict[fastq_name] = {}

                    if "R1" in fastq:
                        fastq_dict[fastq_name]['R1'] = fastq_path
                    if "R2" in fastq:
                        fastq_dict[fastq_name]['R2'] = fastq_path
                        
        #Print fastq_dict
        print(json.dumps(fastq_dict, indent=4))
    
    return fastq_dict

def obtain_final_RBs(file_Rbs, final_RBs):
    

     # Open Moscat-RBs.xlsx to .csv
    df = pd.read_excel(file_Rbs)
    df.to_csv ('Moscat-RBs.csv',
                index= False,
                header = True)
    df_Moscat= pd.DataFrame(pd.read_csv('Moscat-RBs.csv'))
    
    # Obtain all RBs as a list
    RBs_Moscat = df_Moscat['RB code'].tolist()

    #Count how many RBs in 'Moscat-RBs.xlsx' --> 909
    num_RBs_Moscat = len(RBs_Moscat)

    # Compare RBs in Moscat with fastqs names
    existent_RB = []
    non_existent_RB = []
    num_existent_RB = 0
    num_non_existent_RB = 0

    for RB in RBs_Moscat:

        if RB in fastq_dict:
            existent_RB.extend([RB])
            print(f'{RB} fastq exists and its path is :{fastq_dict[RB]}º\n')
            num_existent_RB += 1       

        else:
            print(f'{RB} fastq does not exist') 
            non_existent_RB.extend([RB])
            num_non_existent_RB +=1
        
    
    print(num_RBs_Moscat)       
    print(num_existent_RB)
    print(num_non_existent_RB)

    final_RBs = []
    num_final_RBs = 0
    for RBs in existent_RB:
        fastq_R1_R2 = fastq_dict.get(RBs)
        if 'R1' in fastq_R1_R2 and 'R2' in fastq_R1_R2:
            final_RBs.extend([RBs])
            num_final_RBs +=1

    print(num_final_RBs)
    print(final_RBs)

    return final_RBs

def download_fastqs(fastq_dir, final_RBs, fastq_dict) :
        
    download_log = os.path.join(fastq_dir, "downloads.log")
    o = open(download_log, "w")
    with ftputil.FTPHost("x", "x", "x") as ftp_host:
        for RB in final_RBs:
            paths = fastq_dict.get(RB)

            #Obtain R1 and R2 paths 
            path_R1 = paths.get('R1')
            path_R2 = paths.get('R2')

        
            #Obtain the name of the files for R1 and R2
            fq1_name = os.path.basename(path_R1)
            fq2_name = os.path.basename(path_R2)

            # output fastq
            output_fq1 = os.path.join(fastq_dir, fq1_name)
            output_fq2 = os.path.join(fastq_dir, fq2_name)

            # Download fq1
            if ftp_host.path.isfile(path_R1):
                msg = f" INFO: Descarregant fq1 {path_R1} a {output_fq1}"
                print(msg)
                o.write(msg+"\n")
                ftp_host.download(path_R1, output_fq1)

            # Download fq2
            if ftp_host.path.isfile(path_R2):
                msg = f" INFO: Descarregant fq2 {path_R2} a {output_fq2}"
                print(msg)
                o.write(msg+"\n")
                ftp_host.download(path_R2, output_fq2)
                 
    o.close()







if __name__ == "__main__":

    synology_path = "/CGC/FASTQ_Files_MiSeq/CLINICA/" 
    fastq_dir = "/home/gencardio/tfm_eli/all_files"
    file_Rbs = 'Moscat-RBs.xlsx' 

    fastq_dict= dict()
    fastq_dict = fastqs_identification(synology_path, fastq_dict)

    final_RBs = []
    final_RBs = obtain_final_RBs(file_Rbs,final_RBs)   

    #download_fastqs(fastq_dir, final_RBs, fastq_dict)  


