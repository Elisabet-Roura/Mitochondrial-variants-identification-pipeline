#!/usr/bin/env python

import sys
import os
import subprocess
import logging

#Function to perform variant calling for chrM and chrM shifted files

def GATK_module(removed_duplicates_files_list, gatk_version, genome_reference, genome_reference_shifted, reference_path, reference_path_shifted):

    vcf_files_list = []


    for removed_duplicates_file in removed_duplicates_files_list:

        variant_calling_path_original= os.path.dirname(removed_duplicates_file).split('chrM')[0]
        variant_calling_path = variant_calling_path_original.replace('mapping','variant_calling') + genome_reference
        variant_calling_path_shifted = variant_calling_path_original.replace('mapping','variant_calling') + genome_reference_shifted


        if not os.path.exists(variant_calling_path):
            os.makedirs(variant_calling_path)
        if not os.path.exists(variant_calling_path_shifted):
            os.makedirs(variant_calling_path_shifted)
            
    
        BAM_FOLDER = os.path.dirname(removed_duplicates_file)
        READY_BAM_NAME = os.path.basename(removed_duplicates_file)

        VCF_FOLDER = variant_calling_path
        VCF_FOLDER_SHIFTED = variant_calling_path_shifted
        
        MUTECT2_VCF_NAME = READY_BAM_NAME.replace('.rd.bam', '.mutect2.vcf.gz')

        REF_DIRNAME = os.path.dirname(reference_path)
        REF_FASTA_NAME = os.path.basename(reference_path)

        REF_DIRNAME_SHIFTED = os.path.dirname(reference_path_shifted)
        REF_FASTA_NAME_SHIFTED = os.path.basename(reference_path_shifted)

        sample_name = READY_BAM_NAME.split('_')[0]

        vcf_final_path = os.path.join(VCF_FOLDER,MUTECT2_VCF_NAME)
        vcf_final_path_shifted = os.path.join(VCF_FOLDER_SHIFTED,MUTECT2_VCF_NAME)

        #chrM

        if 'shifted' not in READY_BAM_NAME:

            if os.path.exists(vcf_final_path):
                vcf_files_list.append(vcf_final_path)

            if not os.path.exists(vcf_final_path):

                cmd_mutect2 = f'docker run -v {BAM_FOLDER}:/bam_data/\
                -v {VCF_FOLDER}:/vcf_data/\
                -v {REF_DIRNAME}:/bundle/\
                -it broadinstitute/gatk:{gatk_version} \
                gatk Mutect2 -R /bundle/{REF_FASTA_NAME}\
                -L chrM   --mitochondria-mode\
                -I /bam_data/{READY_BAM_NAME}\
                -O /vcf_data/{MUTECT2_VCF_NAME}'
                p1 = subprocess.run(cmd_mutect2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(cmd_mutect2)
                

                output = p1.stdout.decode("UTF-8")
                error = p1.stderr.decode("UTF-8")
                if p1.returncode != 0:
                    msg = f" ERROR: Could not run Mutect2 for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error)
                    sys.exit()
                
                if p1.returncode == 0:
                    vcf_files_list.append(vcf_final_path)
        

        # chrM shifted

        else:

            if os.path.exists(vcf_final_path_shifted):
                vcf_files_list.append(vcf_final_path_shifted)

            if not os.path.exists(vcf_final_path_shifted):

                cmd_mutect2 = f'docker run -v {BAM_FOLDER}:/bam_data/\
                -v {VCF_FOLDER_SHIFTED}:/vcf_data/\
                -v {REF_DIRNAME_SHIFTED}:/bundle/\
                -it broadinstitute/gatk:{gatk_version} \
                gatk Mutect2 -R /bundle/{REF_FASTA_NAME_SHIFTED}\
                -L chrM   --mitochondria-mode\
                -I /bam_data/{READY_BAM_NAME}\
                -O /vcf_data/{MUTECT2_VCF_NAME}'
        
                p2 = subprocess.run(cmd_mutect2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(cmd_mutect2)
                

                output = p2.stdout.decode("UTF-8")
                error = p2.stderr.decode("UTF-8")
                if p2.returncode != 0:
                    msg = f" ERROR: Could not run Mutect2 for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error)
                    sys.exit()
                
                if p2.returncode == 0:
                    vcf_files_list.append(vcf_final_path_shifted)
            
        
        
    print(vcf_files_list)
    
    return vcf_files_list

#Function to Liftover the output VCF files

def liftovervcf(vcf_files_list, picard_path,reference_path, chain_file):

    lift_over_vcf_list = []


    for vcf_file in vcf_files_list:

        if 'shifted' in vcf_file:

            sample_name = os.path.basename(vcf_file).split('_')[0]
            fixed_vcf = vcf_file.replace('.mutect2.vcf.gz','.fixed.mutect2.vcf.gz')
            rejected_variants = vcf_file.replace('.mutect2.vcf.gz','_rejected_variants.vcf.gz')

            if os.path.exists(fixed_vcf):

                lift_over_vcf_list.append(fixed_vcf)
            
            else:

                cmd_liftover = f'java -jar {picard_path} LiftoverVcf\
                -I {vcf_file} \
                -O {fixed_vcf}\
                -C {chain_file} \
                -R {reference_path}\
                --REJECT {rejected_variants}'

                p1 = subprocess.run(cmd_liftover, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(cmd_liftover)

                output = p1.stdout.decode("UTF-8")
                error = p1.stderr.decode("UTF-8")
                if p1.returncode != 0:
                    msg = f" ERROR: Could not run LiftOverVCFfor sample {sample_name}"
                    logging.error(msg)
                    logging.error(error)
                    sys.exit()
                
                if p1.returncode == 0:
                    lift_over_vcf_list.append(fixed_vcf)

        else:
            lift_over_vcf_list.append(vcf_file)
    

    print(lift_over_vcf_list)
    return(lift_over_vcf_list)


#Function to combine the variant calls from the control region with the non-control region

def merge_vcf (lift_over_vcf_list, picard_path):

    merged_vcf_list = []

    for lift_over_vcf_file in lift_over_vcf_list:

        if not 'shifted' in lift_over_vcf_file:
                
            sample_name = os.path.basename(lift_over_vcf_file).split('_')[0]
            final_merged_path = os.path.dirname(lift_over_vcf_file).split('/chrM')[0]
            final_merged_file_name = '/' + sample_name + '_merged.vcf.gz'
            final_merged_vcf = final_merged_path + final_merged_file_name
            lift_over_file_2 = lift_over_vcf_file.replace('chrM','chrM_shifted').replace('.mutect2.vcf.gz','.fixed.mutect2.vcf.gz')
    

            if os.path.exists(final_merged_vcf):
                merged_vcf_list.append(final_merged_vcf)

            else: 

                cmd_merge_vcfs = f'java -jar {picard_path} MergeVcfs\
                -I {lift_over_vcf_file} \
                -I {lift_over_file_2}\
                -O {final_merged_vcf}'
        
                p1 = subprocess.run(cmd_merge_vcfs, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(cmd_merge_vcfs)

                output = p1.stdout.decode("UTF-8")
                error = p1.stderr.decode("UTF-8")
                if p1.returncode != 0:
                    msg = f" ERROR: Could not run Merge Vcfs for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error)
                    sys.exit()
                
                if p1.returncode == 0:
                    merged_vcf_list.append(final_merged_vcf)
 
    print(merged_vcf_list)
    return merged_vcf_list

# Function to merge stats files for output VCFs

def merge_vcfs_stats(vcf_files_list, gatk_version):

    merged_stats_list = []

    for vcf_file in vcf_files_list:

        if not 'shifted' in vcf_file:

            vcf_stats_path = os.path.dirname(vcf_file)
            vcf_stats_path_shifted = vcf_stats_path.replace('chrM', 'chrM_shifted')
            vcf_stats_path_merged = vcf_stats_path.split('/chrM')[0]

            sample_name = os.path.basename(vcf_file).split('.')[0]

            vcf_stats_file = os.path.basename(vcf_file).replace('.vcf.gz','.vcf.gz.stats')
            vcf_stats_file_shifted = vcf_stats_file.replace('chrM', 'chrM_shifted')
            vcf_stats_file_merged = sample_name + '.vcf.gz.stats'

            vcf_stats_merged = vcf_stats_path_merged + '/' + vcf_stats_file_merged 


            if not os.path.exists(vcf_stats_merged):

                cmd_merge_stats = f'docker run -v {vcf_stats_path}:/stats_data/\
                -v {vcf_stats_path_shifted}:/stats_shifted_data/\
                -v {vcf_stats_path_merged}:/stats_merged/\
                -it broadinstitute/gatk:{gatk_version} \
                gatk MergeMutectStats -stats /stats_data/{vcf_stats_file} \
                -stats /stats_shifted_data/{vcf_stats_file_shifted}\
                -O /stats_merged/{vcf_stats_file_merged}'
        
                p1 = subprocess.run(cmd_merge_stats, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(cmd_merge_stats)

                output = p1.stdout.decode("UTF-8")
                error = p1.stderr.decode("UTF-8")
                if p1.returncode != 0:
                    msg = f" ERROR: Could not run Merge stats for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error)
                    sys.exit()
                
                if p1.returncode == 0 and not vcf_stats_merged in merged_stats_list:
                    merged_stats_list.append(vcf_stats_merged)
            
            else:
                if os.path.exists(vcf_stats_merged) and not vcf_stats_merged in merged_stats_list:
                    merged_stats_list.append(vcf_stats_merged)
    
    print(merged_stats_list)
    return(merged_stats_list)


#FUnction to filter the Variant Calls by Parameters

def filter_mutect_calls(vcf_files_list, gatk_version, reference_path):

    vcf_filtered_list = []


    for vcf_file in vcf_files_list: 
        

        vcf_file_path = os.path.dirname(vcf_file)
        vcf_file_name = os.path.basename(vcf_file)
        sample_name = vcf_file_name.split('.')[0]
        vcf_file_filtered = vcf_file_name.replace('.vcf.gz','_filtered.vcf.gz')
        final_filtered_path = os.path.join(vcf_file_path,vcf_file_filtered)


        reference_dir = os.path.dirname(reference_path)
        reference_file = os.path.basename(reference_path)

        stats_file = vcf_file_name.replace('.chrM_merged.vcf.gz','.vcf.gz.stats')

        if not os.path.exists(final_filtered_path):

            cmd_filter = f'docker run -v {vcf_file_path}:/vcf_data/\
            -v {reference_dir}:/reference_data/\
            -it broadinstitute/gatk:{gatk_version} \
            gatk FilterMutectCalls -R /reference_data/{reference_file} \
            -V /vcf_data/{vcf_file_name}\
            --stats /vcf_data/{stats_file}\
            -O /vcf_data/{vcf_file_filtered}\
            --mitochondria-mode'
            # --autosomal-coverage\ -- NO 

            p1 = subprocess.run(cmd_filter, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(cmd_filter)

            output = p1.stdout.decode("UTF-8")
            error = p1.stderr.decode("UTF-8")
            if p1.returncode != 0:
                msg = f" ERROR: Could not run Filter Mutect Calls for sample {sample_name}"
                logging.error(msg)
                logging.error(error)
                sys.exit()
            
            if p1.returncode == 0 and not final_filtered_path in vcf_filtered_list:
                vcf_filtered_list.append(final_filtered_path)
        
        else:
            if not final_filtered_path in vcf_filtered_list:
                vcf_filtered_list.append(final_filtered_path)
    
    print(vcf_filtered_list)
    return(vcf_filtered_list)



#As a module

if __name__ == "__main__":
    GATK_module()
    liftovervcf()
    merge_vcf()
    merge_vcfs_stats()
    filter_mutect_calls()