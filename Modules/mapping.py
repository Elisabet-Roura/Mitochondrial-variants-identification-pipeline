#!/usr/bin/env python

import sys
import os
import subprocess
import logging



def download_reference_genome(download_genome_dir,genome_version):
    
    #Create the directory to downlad the genome

    if not os.path.isdir(download_genome_dir):
        os.mkdir(download_genome_dir)
    
    os.chdir(download_genome_dir)
    
    if not os.path.exists(f'{download_genome_dir}/{genome_version}.fa.gz'):

        #Download the genome
        download = f'wget https://hgdownload.soe.ucsc.edu/goldenPath/{genome_version}/bigZips/{genome_version}.fa.gz'

        #Unzip de genome
        unzip = f'gunzip {genome_version}.fa.gz'

        #Index the genome
        index = f'bwa index {genome_version}.fa'

        subprocess.run(download, shell = True)
        subprocess.run(unzip, shell = True)
        subprocess.run(index, shell = True)
        print(download, unzip, index)
    
    genome_path = f'{download_genome_dir}/{genome_version}.fa'

    
    return genome_path


def download_reference_chromosome(download_genome_dir,genome_version,chromosome_name):
    
    #Create the directory to downlad the genome:
     
    if not os.path.isdir(download_genome_dir):
        os.mkdir(download_genome_dir)
    
    os.chdir(download_genome_dir)
    
    if not os.path.exists(f'{download_genome_dir}/{chromosome_name}.fa.gz'):

        #Download the genome
        download = f'wget https://hgdownload.soe.ucsc.edu/goldenPath/{genome_version}/chromosomes/{chromosome_name}.fa.gz'

        #Unzip de genome
        unzip = f'gunzip {chromosome_name}.fa.gz'

        #Index the genome
        index = f'bwa index {chromosome_name}.fa'

        subprocess.run(download, shell = True)
        subprocess.run(unzip, shell = True)
        subprocess.run(index, shell = True)
        print(download, unzip, index)
    
    chrM_path = f'{download_genome_dir}/{chromosome_name}.fa'

    
    return chrM_path


#Function to obtain the sorted.bam files name

def sorted_bam_files(input_fastq, genome_reference):

    final_sorted_bam_list = []

    fastqs_dir = os.listdir(input_fastq)
    print(fastqs_dir)

    for sample_name in fastqs_dir:
        if sample_name.startswith("RB"):

            sample_dir = os.path.join(input_fastq, sample_name)
            sample_dir = sample_dir +'/mapping'

            final_file_directory = f'{sample_dir}/{sample_name}_{genome_reference}_sorted.bam'
            final_sorted_bam_list.append(final_file_directory)
    
    print(final_sorted_bam_list)

    return final_sorted_bam_list


#Function to map the fastqs to obtain .sam file

def mapping_fastqs(input_fastq, paired_dictionary, genome_path, genome_reference):


    sam_files_list = []
    print(genome_path)
    
    fastqs_dir = os.listdir(input_fastq)
    
    for sample_name in fastqs_dir:

        if sample_name.startswith("RB"):

            sample_dir = os.path.join(input_fastq, sample_name)
            sample_dir = sample_dir +'/mapping/'+ genome_reference

            if not os.path.isdir(sample_dir):
                os.makedirs(sample_dir)
        
        sam_path = f'{sample_dir}/{sample_name}_{genome_reference}.sam'
        sorted_bam_path = f'{sample_dir}/{sample_name}_{genome_reference}_sorted.bam'
        
        if not os.path.exists(sorted_bam_path):     

            if sample_name in paired_dictionary:
            
                
                trimmed_R1_file = paired_dictionary[sample_name]['TRIMMED_R1'] 
                trimmed_R2_file = paired_dictionary[sample_name]['TRIMMED_R2'] 

                print(trimmed_R1_file,trimmed_R2_file)
                print(genome_path)

                    
                if not os.path.exists(sam_path):

                    cmd = f'bwa mem {genome_path} {trimmed_R1_file} {trimmed_R2_file} -t 4 > {sam_path}'
                    print(cmd)
                    p_map = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    
                    output_map = p_map.stdout.decode("UTF-8")
                    error_map = p_map.stderr.decode("UTF-8")

                    if p_map.returncode != 0:
                        msg = f"ERROR: Could not run mapping command for sample {sample_name}"
                        logging.error(msg)
                        logging.error(error_map)
                        sys.exit(1)

                    else:
                        sam_files_list.append(sam_path)

                else:
                    if os.path.exists(sam_path) and not sam_path in sam_files_list:
                            sam_files_list.append(sam_path) 
        
        if os.path.exists(sorted_bam_path) and not sam_path in sam_files_list:
            sam_files_list.append(sam_path)


    print(sam_files_list)
         
    return sam_files_list


#Funtcion to obtain a unsorted.bam form each .sam

def from_sam_to_bam(sam_files_list, final_sorted_bam_files):

    unsorted_bam_files_list = []
    
    for sorted_bam_file in final_sorted_bam_files:

        sam_path = sorted_bam_file.replace("_sorted.bam", ".sam")

        if sam_path in sam_files_list:

            unsorted_bam_path =  sam_path.replace(".sam","_unsorted.bam")
            sample_name = os.path.basename(unsorted_bam_path).split('_')[0]


            if os.path.exists(sorted_bam_file):
                if not unsorted_bam_path in unsorted_bam_files_list:
                    unsorted_bam_files_list.append(unsorted_bam_path)

            if not os.path.exists(sorted_bam_file):
  
                if os.path.exists(sam_path):
                    
                    if not os.path.exists(unsorted_bam_path):
                        if not unsorted_bam_path in unsorted_bam_files_list: 
                            
                            cmd = f'samtools view -bS {sam_path} > {unsorted_bam_path}'
                            print(cmd)
                            p_bam = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            

                            output_bam = p_bam.stdout.decode("UTF-8")
                            error_bam = p_bam.stderr.decode("UTF-8")

                            if p_bam.returncode != 0:
                                msg = f"ERROR: Could not convert SAM to BAM for sample {sample_name}"
                                logging.error(msg)
                                logging.error(error_bam)
                                sys.exit(1)
                            else:
                                unsorted_bam_files_list.append(unsorted_bam_path)

            
                if os.path.exists(unsorted_bam_path):
                    if not unsorted_bam_path in unsorted_bam_files_list:
                        unsorted_bam_files_list.append(unsorted_bam_path)



                        
    print(unsorted_bam_files_list)

    return unsorted_bam_files_list
        

#Function to sort de unsorted.bam to obtain a sorted.bam

def from_unsorted_to_sorted_bam(unsorted_bam_files_list):

    sorted_bam_files_list = []

    for unsorted_bam_file in unsorted_bam_files_list:

        sorted_bam_file = unsorted_bam_file.replace("_unsorted.bam", "_sorted.bam")
        sample_name = os.path.basename(unsorted_bam_file).split('_')[0]


        if not os.path.exists(sorted_bam_file): 

            cmd = f'samtools sort -T TEST {unsorted_bam_file} -o {sorted_bam_file}'
            print(cmd)
            p_sort = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output_sort = p_sort.stdout.decode("UTF-8")
            error_sort = p_sort.stderr.decode("UTF-8")

            if p_sort.returncode != 0:
                msg = f"ERROR: Could not sort BAM file for sample {sample_name}"
                logging.error(msg)
                logging.error(error_sort)
                sys.exit(1)
            else:
                sorted_bam_files_list.append(sorted_bam_file)

        if os.path.exists(sorted_bam_file): 
            if not sorted_bam_file in sorted_bam_files_list:
                sorted_bam_files_list.append(sorted_bam_file)
    
        
    
    print(sorted_bam_files_list)

    return sorted_bam_files_list

#Function to eliminate .sam and unsorted.bam files

        
def remove_files(sorted_bam_files_list):

    for sorted_bam_file in sorted_bam_files_list:
        if os.path.exists(sorted_bam_file):

            #Eliminate .sam
            sam_file_dir = sorted_bam_file.replace('_sorted.bam','.sam')

            if not os.path.exists(sam_file_dir):
                print(f'{sam_file_dir} does not exist')
            
            else: 
                if os.path.exists(sam_file_dir):
                    os.remove(sam_file_dir)
                    msg = f'{sam_file_dir} has been deleted'
                    print(msg)

            
            
            #Eliminate .unsorted.bam
            unsorted_bam_dir = sorted_bam_file.replace('_sorted.bam','_unsorted.bam')

            if not os.path.exists(unsorted_bam_dir):
                print(f'{unsorted_bam_dir} does not exist')

            else:
                if os.path.exists(unsorted_bam_dir):
                    os.remove(unsorted_bam_dir)
                    msg = f'{unsorted_bam_dir} has been deleted'
                    print(msg)

        
# mapp_sort_eliminate function combines mapping_fastqs, from_sam_to_bam, from_unsorted_to_sorted_bam, and remove_files; 
# you can either use this single function for the mapping step or run the smaller functions sequentially.


def mapp_sort_eliminate(final_sorted_files_list, paired_dictionary,  genome_path):

    sorted_bam_files_list = []

    for final_sorted_file in final_sorted_files_list:

        sorted_file_path = os.path.dirname(final_sorted_file)
        index_bam = final_sorted_file.replace('.bam','bam.bai')

        if not os.path.isdir(sorted_file_path):
            os.makedirs(sorted_file_path)
        
        if not os.path.exists(final_sorted_file):

            sam_path = final_sorted_file.replace('_sorted.bam','.sam')
            unsorted_bam_path = final_sorted_file.replace('_sorted.bam', '_unsorted.bam')
            sample_name = os.path.basename(final_sorted_file).split('_')[0]

            if sample_name in paired_dictionary:
                
                trimmed_R1_file = paired_dictionary[sample_name]['TRIMMED_R1'] 
                trimmed_R2_file = paired_dictionary[sample_name]['TRIMMED_R2'] 

                print(trimmed_R1_file,trimmed_R2_file)
                print(genome_path)

                #Map fastqs
                    
                if not os.path.exists(sam_path):

                    cmd = f'bwa mem {genome_path} {trimmed_R1_file} {trimmed_R2_file} -t 4 > {sam_path}'
                    print(cmd)
                    p_map = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    
                    output_map = p_map.stdout.decode("UTF-8")
                    error_map = p_map.stderr.decode("UTF-8")

                    if p_map.returncode != 0:
                        msg = f"ERROR: Could not run mapping command for sample {sample_name}"
                        logging.error(msg)
                        logging.error(error_map)
                        sys.exit(1)
                
                #From .sam to _unsorted.bam

                if not os.path.exists(unsorted_bam_path):
                    
                    cmd = f'samtools view -bS {sam_path} > {unsorted_bam_path}'
                    print(cmd)
                    p_bam = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    

                    output_bam = p_bam.stdout.decode("UTF-8")
                    error_bam = p_bam.stderr.decode("UTF-8")

                    if p_bam.returncode != 0:
                        msg = f"ERROR: Could not convert SAM to BAM for sample {sample_name}"
                        logging.error(msg)
                        logging.error(error_bam)
                        sys.exit(1)

                #From _unsorted.bam to _sorted.bam

                cmd = f'samtools sort -T TEST {unsorted_bam_path} -o {final_sorted_file}'
                cmd_2 = f' samtools index {final_sorted_file}'
                print(cmd, cmd_2)
                p_sort = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p_index = subprocess.run(cmd_2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                

                output_sort = p_sort.stdout.decode("UTF-8")
                error_sort = p_sort.stderr.decode("UTF-8")

                if p_sort.returncode != 0:
                    msg = f"ERROR: Could not sort BAM file for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_sort)
                    sys.exit(1)
                else:
                    sorted_bam_files_list.append(final_sorted_file)
            else:
                print(f'{sample_name} not in dictionary')
        

        #Cleanup unecessary files after mapping 

        if os.path.exists(final_sorted_file):

            if not os.path.exists(index_bam):
                cmd_2 = f' samtools index {final_sorted_file}'
                print(cmd_2)
                p_index = subprocess.run(cmd_2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)



            #Eliminate .sam
            sam_file_dir = final_sorted_file.replace('_sorted.bam','.sam')

            if not os.path.exists(sam_file_dir):
                print(f'{sam_file_dir} does not exist')
            
            else: 
                os.remove(sam_file_dir)
                msg = f'{sam_file_dir} has been deleted'
                print(msg)

            
            #Eliminate .unsorted.bam
            unsorted_bam_dir = final_sorted_file.replace('_sorted.bam','_unsorted.bam')

            if not os.path.exists(unsorted_bam_dir):
                print(f'{unsorted_bam_dir} does not exist')

            else:

                os.remove(unsorted_bam_dir)
                msg = f'{unsorted_bam_dir} has been deleted'
                print(msg)
            

            if not final_sorted_file in sorted_bam_files_list:
                
                sorted_bam_files_list.append(final_sorted_file)            
    


    print(sorted_bam_files_list)  
    return(sorted_bam_files_list)

    
    
#Function to preprocess sorted.bam files before variant_calling calling extract_chrM.pl

def preprocessing_variant_calling (sorted_bam_files_list, perl_script):

    merged_bam_alignment_files_list = []

    for bam in sorted_bam_files_list:

        bam_dir = os.path.dirname(bam)
        bam_file = os.path.basename(bam)
        sam_path = bam.replace('_hg38_sorted.bam','_chrM.sam')
        merged_bam_name = bam_file.replace('_hg38_sorted.bam', '.chrM_merged.bam')
        merged_bam_name_shifted= bam_file.replace('_hg38_sorted.bam', '.chrM_shifted_merged.bam')
        merged_bam_file = bam_dir +'/chrM/'+ merged_bam_name
        merged_bam_shifted_file = bam_dir +'/chrM_shifted/'+ merged_bam_name_shifted


        if not os.path.exists(merged_bam_file) and not os.path.exists(merged_bam_shifted_file):
            cmd = ["perl", perl_script, bam]
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            print(f"Ejecutado para: {bam}")
            print(result.stdout.decode())

            if result.returncode != 0:
                print(f"ERROR para {bam}:\n{result.stderr.decode()}")
            
            else:
                if not merged_bam_file in merged_bam_alignment_files_list:
                    merged_bam_alignment_files_list.append(merged_bam_file)
                if not merged_bam_shifted_file in merged_bam_alignment_files_list:
                    merged_bam_alignment_files_list.append(merged_bam_shifted_file)
                
                if os.path.exists(sam_path):
                    os.remove(sam_path)
                    print(f'{sam_path} has been deleted')

        if os.path.exists(merged_bam_file) and os.path.exists(merged_bam_shifted_file): 
            if not merged_bam_file in merged_bam_alignment_files_list:
                merged_bam_alignment_files_list.append(merged_bam_file)
            if not merged_bam_shifted_file in merged_bam_alignment_files_list:
                merged_bam_alignment_files_list.append(merged_bam_shifted_file)

            if os.path.exists(sam_path):
                os.remove(sam_path)
                print(f'{sam_path} has been deleted')

    print(merged_bam_alignment_files_list)
    return(merged_bam_alignment_files_list)


          
#Function to eliminate duplicates from the .bam files, adding before the read groups

def eliminate_duplicates(bam_files_list, picard_path):

    removed_duplicate_files_list = []

    for bam_file in bam_files_list:

        bam_file_name = os.path.basename(bam_file)
        sample_name = bam_file_name.split('_')[0]
        removed_duplicates_file = bam_file.replace('.bam','.rd.bam')
        read_groups_file = bam_file.replace('.bam','.rg.bam')
        metrics_file = bam_file.replace('.bam','.metrics.txt')

        
        if not os.path.exists(removed_duplicates_file):
            if os.path.exists(bam_file):

                cmd_1 = f'java -jar {picard_path}  AddOrReplaceReadGroups -I {bam_file}\
                -O {read_groups_file}\
                -RGID {sample_name}\
                -LB lib1\
                -PL illumina\
                -PU unit1\
                -SM {sample_name} --CREATE_INDEX'
                cmd_2 = f'java -jar {picard_path} MarkDuplicates\
                -I {read_groups_file}\
                -O {removed_duplicates_file}\
                -M {metrics_file} --REMOVE_DUPLICATES True --CREATE_INDEX'

                print(cmd_1)
                p_rg = subprocess.run(cmd_1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                

                output_rg = p_rg.stdout.decode("UTF-8")
                error_rg = p_rg.stderr.decode("UTF-8")

                if p_rg.returncode != 0:
                    msg = f"ERROR: Could not run read group for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_rg)
                    sys.exit(1)
                
               
                    
                print(cmd_2)
                p_md = subprocess.run(cmd_2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                

                output_md = p_md.stdout.decode("UTF-8")
                error_md = p_md.stderr.decode("UTF-8")

                if p_md.returncode != 0:
                    msg = f"ERROR: Could not run mark duplicates for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_md)
                    sys.exit(1)
           
                else:
                    if removed_duplicates_file not in removed_duplicate_files_list:
                        removed_duplicate_files_list.append(removed_duplicates_file)

            if not os.path.exists(bam_file):
                print(f'Failed to remove duplicates for sample {sample_name}')

        if os.path.exists(removed_duplicates_file) and not removed_duplicates_file in removed_duplicate_files_list: 
            removed_duplicate_files_list.append(removed_duplicates_file)
        
    print(removed_duplicate_files_list)   

    return removed_duplicate_files_list



#Function to index .bam files

def index_bam(removed_duplicate_files_list):

    index_bam_files_list = []

    for removed_duplicates_file in removed_duplicate_files_list:

        index_sorted_bam_file = removed_duplicates_file.replace('.bam', '.bam.bai')
        sample_name = os.path.basename(removed_duplicates_file).split('_')[0]

        if not os.path.exists(index_sorted_bam_file):

            cmd = f'samtools index {removed_duplicates_file}'
            print(cmd)
            p_index = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output_index = p_index.stdout.decode("UTF-8")
            error_index = p_index.stderr.decode("UTF-8")

            if p_index.returncode != 0:
                msg = f"ERROR: Could not index removed_duplicates BAM file for sample {sample_name}"
                logging.error(msg)
                logging.error(error_index)
                sys.exit(1)
            else:
                index_bam_files_list.append(index_sorted_bam_file)
                
        if os.path.exists(index_sorted_bam_file):
            if not index_sorted_bam_file in index_bam_files_list:
                index_bam_files_list.append(index_sorted_bam_file)
    
    print(index_bam_files_list)

    return index_bam_files_list   




                  


#As a module

if __name__ == "__main__":

    download_reference_genome()
    download_reference_chromosome()
    mapping_fastqs()
    sorted_bam_files()
    from_sam_to_bam()
    from_unsorted_to_sorted_bam()
    mapp_sort_eliminate()
    preprocessing_variant_calling()
    eliminate_duplicates()
    index_bam()
    remove_files()
  

