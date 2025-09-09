#!/usr/bin/env python

import os
import subprocess


#Function to anotate the vcf_files calling the annotation perl script

def annotation(annotation_script, vcf_final_files):

    annotated_files_list = []

    for final_vcf in vcf_final_files:

        sample_name = os.path.basename(final_vcf).split('.')[0]
        variant_calling_path = os.path.dirname(final_vcf)
        annotation_path = variant_calling_path.replace('variant_calling', 'annotation')
        annotated_file = annotation_path + '/'+ sample_name + '_annotated.vcf'

        if not os.path.exists(annotation_path):
            os.makedirs(annotation_path)

        
        if os.path.exists(annotated_file):
            if not annotated_file in annotated_files_list:
                annotated_files_list.append(annotated_file)
        
        else:

            cmd = [ "perl", annotation_script, final_vcf ]
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            print(f"Ejecutado para: {sample_name}")
            print(result.stdout.decode())

            if result.returncode != 0:
                print(f"ERROR para {sample_name}:\n{result.stderr.decode()}")
            
            else:
                if not annotated_file in annotated_files_list:
                    annotated_files_list.append(annotated_file)

    print(annotated_files_list)
    print(len(annotated_files_list))
    return(annotated_files_list)


if __name__ == "__main__":

    annotation()


        

        
