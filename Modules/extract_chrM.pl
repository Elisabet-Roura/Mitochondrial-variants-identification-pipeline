
#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);

my $input_bam = $ARGV[0];


my $ref_fasta_chrM = "/media/gencardio/KINGSTON/hg38/hg38.fa";
my $ref_fasta_chrM_shifted = "/media/gencardio/KINGSTON/hg38/hg38_with_shifted_chrM.fa";
my $picard  = "/home/gencardio/tfm_eli/PIPELINE/picard.jar";


my $bam_dir =  dirname($input_bam);


my $sample_name = basename($input_bam);
$sample_name =~ s/_.*\.bam$//;



# Input and output files for chrM

my $sam_file  = $bam_dir. "/" .$sample_name . "_chrM.sam";
my $bam_file  =  $bam_dir. "/" .$sample_name .  "_chrM.bam";
my $unaln_bam  =  $bam_dir. "/" .$sample_name . ".unaligned.bam";
my $sort_aln   =  $bam_dir. "/" .$sample_name . ".querysorted.bam";
my $sort_unaln  =  $bam_dir. "/" .$sample_name . ".unaligned.querysorted.bam";
my $merged_bam   =  $bam_dir."/chrM/".$sample_name . ".chrM_merged.bam";

# Input and output files chrM shifted

my $sam_file_shifted  =  $bam_dir. "/" .$sample_name . "_chrM_shifted.sam";
my $bam_file_shifted   =  $bam_dir. "/" .$sample_name .  "_chrM_shifted.bam";
my $unaln_bam_shifted   =  $bam_dir. "/" .$sample_name . ".unaligned_shifted.bam";
my $sort_aln_shifted    =  $bam_dir. "/" .$sample_name . ".querysorted_shifted.bam";
my $sort_unaln_shifted   =  $bam_dir. "/" .$sample_name . ".unaligned.querysorted_shifted.bam";
my $merged_bam_shifted    =  $bam_dir."/chrM_shifted/".$sample_name . ".chrM_shifted_merged.bam";

make_path("$bam_dir/chrM");
make_path("$bam_dir/chrM_shifted");


# =========
# main
# =========

#chrM 

extract_chrM_reads($input_bam, $sam_file);
sam_to_bam($input_bam, $bam_file);
# sam_to_bam($input_bam,$bam_file);
revert_bam($bam_file, $unaln_bam);
sort_bams($bam_file, $unaln_bam, $sort_aln, $sort_unaln);
merge_bams($ref_fasta_chrM, $picard, $sort_aln, $sort_unaln, $merged_bam);

#chrM shifted

# extract_chrM_reads($input_bam, $sam_file_shifted);
# sam_to_bam($sam_file_shifted, $bam_file_shifted);
# revert_bam($bam_file_shifted, $unaln_bam_shifted);
# sort_bams($bam_file_shifted, $unaln_bam_shifted, $sort_aln_shifted, $sort_unaln_shifted);
merge_bams($ref_fasta_chrM_shifted, $picard, $sort_aln, $sort_unaln, $merged_bam_shifted);

################################

#Subroutine to obtain the reads that we are interested in 

sub extract_chrM_reads {

    my ($bam, $out_sam) = @_;
    print "INFO: Extracting chrM reads\n";


    my %keep;

    #Samtools view without filters to retrieve all reads from chrM, including discordant, multimapped, etc.

    open(IN, "  samtools view $bam |");
    while (my $line=<IN>) {


        chomp $line;
        my @fields = split("\t", $line);
        my $qname = $fields[0];
        my $rname = $fields[2];

        if ($rname && $rname eq "chrM") {
            $keep{$qname} = 1;
        }
    }
    close IN;

    #SAM header

    open(OUT, ">", $out_sam) or die "Can't write to $out_sam: $!";
    open(HDR, " samtools view -H $bam |");
    while (my $line=<HDR>) {
        print OUT $line;
    }
    close HDR;

    open(IN, " samtools view $bam |");
    while (my $line=<IN>) {
        chomp $line;
        my @fields = split("\t", $line);
        my $qname = $fields[0];

        if ($keep{$qname}) {
            print OUT "$line\n";
        }
    }
    close IN;
    close OUT;
}


#Subroutine to convert sam to bam 

sub sam_to_bam {

    my ($sam, $bam) = @_;
    print "INFO: Converting $sam to bam\n";

    my $cmd = "samtools view -Sb $sam > $bam";
    system($cmd);

    my cmd = "samtools index $bam";
    system($cmd);


}

# sub sam_to_bam {
#     my ($bam, $new_bam) = @_;   # recibes el BAM original y el nombre nuevo

#     print "INFO: Renombrando $bam a $new_bam\n";
#     rename $bam, $new_bam or die "No pude renombrar $bam a $new_bam: $!\n";

#     print "INFO: Indexando $new_bam\n";
#     my $cmd = "samtools index $new_bam";
#     system($cmd) == 0 or die "Error ejecutando: $cmd\n";
# }





#Subroutine to obtain an unaligned bam

sub revert_bam {

    my ($bam, $out_unaligned) = @_;
    print "INFO: Running RevertSam\n";

    my $cmd = "java -jar $picard RevertSam INPUT=$bam OUTPUT=$out_unaligned ";
    system($cmd);
}


#Subroutine to perform sorting by queryname

sub sort_bams {

    my ($aligned_bam, $unaligned_bam, $out_aligned, $out_unaligned) = @_;
    print "INFO: Sorting bams by read name\n";

    my $cmd = "java -jar $picard SortSam -I $aligned_bam -O $out_aligned -SO queryname";
    system($cmd);

    $cmd = "java -jar $picard SortSam -I $unaligned_bam -O $out_unaligned -SO queryname";
    system($cmd);
}


#Subroutine to merge aligned and unaligned bam

sub merge_bams {

    my ($ref, $picard_jar, $aln_bam, $unaln_bam, $out_bam) = @_;
    print "INFO: Running MergeBamAlignment\n";

    my $cmd = "java -jar $picard_jar MergeBamAlignment " .
        "UNMAPPED_BAM=$unaln_bam " .
           "ALIGNED_BAM=$aln_bam " .
            "OUTPUT=$out_bam " .
            "REFERENCE_SEQUENCE=$ref " .
             "CREATE_INDEX=true";
    system($cmd);

}





#--------  HETEROPLASMIA SIMULATION  ADAPTATION----------
#!/usr/bin/env perl
# use strict;
# use warnings;
# use File::Basename;
# use File::Path qw(make_path);

# my $input_bam = $ARGV[0] or die "Uso: $0 <input.bam>\n";

# my $ref_fasta_chrM          = "/media/gencardio/KINGSTON/hg38/hg38.fa";
# my $ref_fasta_chrM_shifted  = "/media/gencardio/KINGSTON/hg38/hg38_with_shifted_chrM.fa";
# my $picard                  = "/home/gencardio/tfm_eli/PIPELINE/picard.jar";

# my $bam_dir = dirname($input_bam);
# my $sample_name = basename($input_bam);
# $sample_name =~ s/_.*\.bam$//;


# my $sam_file        = "$bam_dir/$sample_name\_chrM.sam";
# my $bam_file        = "$bam_dir/$sample_name\_chrM.bam";
# my $bam_rg          = "$bam_dir/$sample_name\_chrM.rg.bam";
# my $unaln_bam       = "$bam_dir/$sample_name.unaligned.bam";
# my $unaln_rg        = "$bam_dir/$sample_name.unaligned.rg.bam";
# my $unaln_unique    = "$bam_dir/$sample_name.unaligned.unique.bam";
# my $unaln_unique_qn = "$bam_dir/$sample_name.unaligned.unique.qn.bam";
# my $sort_aln        = "$bam_dir/$sample_name.querysorted.bam";
# my $sort_aln_fix    = "$bam_dir/$sample_name.querysorted.fixmate.bam";
# my $sort_unaln      = "$bam_dir/$sample_name.unaligned.querysorted.bam";
# my $merged_bam      = "$bam_dir/chrM/$sample_name.chrM_merged.bam";

# # Shifted
# my $merged_bam_shifted = "$bam_dir/chrM_shifted/$sample_name.chrM_shifted_merged.bam";

# make_path("$bam_dir/chrM");
# make_path("$bam_dir/chrM_shifted");

# # =========
# # main
# # =========


# extract_chrM_reads($input_bam, $sam_file);
# sam_to_bam($sam_file, $bam_file);
# add_read_groups($bam_file, $bam_rg, $sample_name);
# revert_bam($bam_rg, $unaln_bam);
# add_read_groups($unaln_bam, $unaln_rg, $sample_name);
# collapse_unmapped_pairs($unaln_rg, $unaln_unique);
# sort_bams($bam_rg, $unaln_unique, $sort_aln, $sort_unaln);
# fix_mate_info($sort_aln, $sort_aln_fix);
# merge_bams($ref_fasta_chrM,$picard, $sort_aln_fix, $sort_unaln, $merged_bam);
# merge_bams($ref_fasta_chrM_shifted, $picard, $sort_aln_fix, $sort_unaln, $merged_bam_shifted);

# exit 0;

# ################################
# # Subrutinas
# ################################


# sub extract_chrM_reads {
#     my ($bam, $out_sam) = @_;
#     print "INFO: Extracting chrM reads (primary only; max 1xR1+1xR2 per QNAME)\n";

#     my %keep;

#     open(my $IN1, "-|", "samtools view $bam") or die $!;
#     while (my $line = <$IN1>) {
#         chomp $line;
#         my @f = split(/\t/, $line);
#         my ($q, $flag, $rname) = @f[0,1,2];
#         next if ($flag & 0x100) || ($flag & 0x800); 
#         if ($rname && $rname eq "chrM") {
#             $keep{$q} = 1;
#         }
#     }
#     close $IN1;


#     open(my $OUT, ">", $out_sam) or die "Can't write $out_sam: $!";
#     open(my $HDR, "-|", "samtools view -H $bam") or die $!;
#     while (my $h = <$HDR>) { print $OUT $h }
#     close $HDR;


#     my (%seenR1, %seenR2);
#     open(my $IN2, "-|", "samtools view $bam") or die $!;
#     while (my $line = <$IN2>) {
#         chomp $line;
#         my @f = split(/\t/, $line);
#         my ($q, $flag) = @f[0,1];
#         next unless $keep{$q};
#         next if ($flag & 0x100) || ($flag & 0x800); 

#         my $is_r1 = int($flag/64)%2;
#         my $is_r2 = int($flag/128)%2;

#         if ($is_r1 && !$seenR1{$q}) { print $OUT $line, "\n"; $seenR1{$q}=1; next; }
#         if ($is_r2 && !$seenR2{$q}) { print $OUT $line, "\n"; $seenR2{$q}=1; next; }

#     }
#     close $IN2;
#     close $OUT;
# }


# sub sam_to_bam {
#     my ($sam, $bam) = @_;
#     print "INFO: Converting $sam -> $bam\n";
#     my $cmd = "samtools view -b -o $bam $sam";
#     system($cmd) == 0 or die "Error: $cmd\n";
#     $cmd = "samtools index $bam";
#     system($cmd) == 0 or die "Error: $cmd\n";
# }

# sub add_read_groups {
#     my ($in_bam, $out_bam, $sm) = @_;
#     print "INFO: Adding Read Groups to $in_bam\n";
#     my $cmd = "java -jar $picard AddOrReplaceReadGroups ".
#               "I=$in_bam O=$out_bam RGID=$sm RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=$sm";
#     system($cmd) == 0 or die "Error: $cmd\n";
# }


# sub revert_bam {
#     my ($bam, $out_unaligned) = @_;
#     print "INFO: Running RevertSam\n";
#     my $cmd = "java -jar $picard RevertSam INPUT=$bam OUTPUT=$out_unaligned";
#     system($cmd) == 0 or die "Error: $cmd\n";
# }

# sub collapse_unmapped_pairs {
#     my ($in_bam, $out_bam) = @_;
#     print "INFO: Collapsing unmapped to 1xR1+1xR2 per QNAME\n";
#     my $header = "$out_bam.header.sam";
#     my $body   = "$out_bam.body.sam";
#     my $cmd;

#     $cmd = "samtools view -H $in_bam > $header";
#     system($cmd) == 0 or die "Error: $cmd\n";

#     $cmd = "samtools view $in_bam | awk -F'\\t' '{q=\$1; f=\$2; r1=int(f/64)%2; r2=int(f/128)%2; ".
#            "if (r1==1){ if(!(q in s1)){ print; s1[q]=1 } } else if (r2==1){ if(!(q in s2)){ print; s2[q]=1 } }}' > $body";
#     system($cmd) == 0 or die "Error: $cmd\n";

#     $cmd = "cat $header $body | samtools view -b -o $out_bam -";
#     system($cmd) == 0 or die "Error: $cmd\n";

#     unlink $header, $body;
# }


# sub sort_bams {
#     my ($aligned_bam, $unaligned_bam, $out_aligned, $out_unaligned) = @_;
#     print "INFO: Sorting bams by read name\n";
#     my $cmd = "java -jar $picard SortSam -I $aligned_bam  -O $out_aligned  -SO queryname";
#     system($cmd) == 0 or die "Error: $cmd\n";
#     $cmd = "java -jar $picard SortSam -I $unaligned_bam -O $out_unaligned -SO queryname";
#     system($cmd) == 0 or die "Error: $cmd\n";
# }


# sub fix_mate_info {
#     my ($in_qn, $out_qn_fix) = @_;
#     print "INFO: Fixing mate info (samtools fixmate)\n";
#     my $cmd = "samtools fixmate -r -m $in_qn $out_qn_fix";
#     system($cmd) == 0 or die "Error: $cmd\n";
# }

# sub merge_bams {
#     my ($ref, $picard_jar, $aln_bam, $unaln_bam, $out_bam) = @_;
#     print "INFO: Running MergeBamAlignment -> $out_bam\n";
#     my $cmd = "java -jar $picard_jar MergeBamAlignment ".
#               "UNMAPPED_BAM=$unaln_bam ".
#               "ALIGNED_BAM=$aln_bam ".
#               "OUTPUT=$out_bam ".
#               "REFERENCE_SEQUENCE=$ref ".
#               "CREATE_INDEX=true PAIRED_RUN=true";
#     system($cmd) == 0 or die "Error: $cmd\n";
# }
