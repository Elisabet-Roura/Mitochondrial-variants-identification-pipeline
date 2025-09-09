#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;




my $input_vcf = $ARGV[0];

if (!$input_vcf) {
    die "Usage: $0 <input_vcf>\n";
}

if(!-e $input_vcf) {
    die "Input VCF file does not exist: $input_vcf\n";
}


# Hypervariable regions file
my $HV = "/home/gencardio/tfm_eli/PIPELINE/modules/Annotation/HV.bed.gz";

# Homoplymeric regions
my $HP = "/home/gencardio/tfm_eli/PIPELINE/modules/Annotation/HP.bed.gz";

# Hotspot regions
my $HS = "/home/gencardio/tfm_eli/PIPELINE/modules/Annotation/HS.bed.gz";

# CDS regions
my $CDS = "/home/gencardio/tfm_eli/PIPELINE/modules/Annotation/CDS.bed.gz";

# Haplogroup database
my $HG = "/home/gencardio/tfm_eli/PIPELINE/modules/Annotation/HG.vcf.gz";

my $gencode_gff = "/home/gencardio/tfm_eli/PIPELINE/modules/Annotation/gencode.chrM.gff.gz";


##############################


# read haplogproups vcf and return a dict
my $haplogroups = vcf2dict($HG);
my %HoHP = %{$haplogroups};


# open input VCF file (From mutect2)

my $annotated_vcf_path = dirname($input_vcf);
$annotated_vcf_path =~s/variant_calling/annotation/;


my $annotated_vcf_name = basename($input_vcf);

$annotated_vcf_name=~ s/\..*//;
$annotated_vcf_name = $annotated_vcf_name . '_annotated.vcf';


my $annotated_vcf = "$annotated_vcf_path/$annotated_vcf_name";


open IN, " zcat $input_vcf | " or die "Cannot open $input_vcf: $!";
open OUT, ">", $annotated_vcf;

while (my $line=<IN>) {
    chomp $line;
    # print "$line\n";
    if ($line=~/#/) {
        print OUT "$line\n";  # Print header lines as is
        next;
    }
    my @tmp = split("\t", $line);
    my $coordinate = "chrM\t$tmp[1]\t$tmp[3]\t$tmp[4]";


    # Haplogroups annotation
    my $haplogroup = ".";
    if (exists $HoHP{$coordinate}) {
        print $HoHP{$coordinate} ."\n";
        $haplogroup = $HoHP{$coordinate};
        $haplogroup =~ s/HG=//;  # Remove the HG prefix
    }

    my $tabixCoordinate = "chrM:$tmp[1]-$tmp[1]";

    # Hotspot regions annotation
    my $hsOut = `tabix $HS $tabixCoordinate`;
    chomp $hsOut;

    my @hsArr = ();

    if ($hsOut) {
        my @hsTmp = split("\n", $hsOut);
        foreach my $ann(@hsTmp) {
            my @tmpAnn = split("\t", $ann);
            my $annValue = $tmpAnn[1] . "_" . $tmpAnn[2];
            push @hsArr, $annValue;
        }
    }
    my $hsAnnotation = join(",", @hsArr);

    # Homopolymeric annotations
    my $hpOut = `tabix $HP $tabixCoordinate`;
    chomp $hpOut;

    my @hpArr = ();
    if ($hpOut) {
        my @hpTmp = split("\n", $hpOut);
        foreach my $ann(@hpTmp) {
            my @tmpAnn = split("\t", $ann);
            my $annValue = $tmpAnn[1] . "_" . $tmpAnn[2];
            push @hpArr, $annValue;
        }
    }
    my $hpAnnotation = join(",", @hpArr);

    # Hypervariable regions annotation
    my $hvOut = `tabix $HV $tabixCoordinate`;
    chomp $hvOut;
    my @hvArr = ();
    if ($hvOut) {
        my @hvTmp = split("\n", $hvOut);
        foreach my $ann(@hvTmp) {
            my @tmpAnn = split("\t", $ann);
            my $annValue = $tmp[3];
            push @hvArr, $annValue;
        }
    }
    my $hvAnnotation = join(",", @hvArr);
    

    # gff annotation
    my $gffOut = `tabix $gencode_gff $tabixCoordinate`;
    chomp $gffOut;


    my @gffArr = ();
    my $gffAnnotation = ".";
    if ($gffOut) {
        my @gffTmp = split("\n", $gffOut);
        foreach my $ann(@gffTmp) {
            my @tmpAnn = split("\t", $ann);

            my $feature = $tmpAnn[2];
            next if $feature ne "exon";

            my @featureInfo = split(";", $tmpAnn[8]);

            my ($gene_name) = grep { /^gene_name=/ } @featureInfo;
            my ($strand) = $tmpAnn[6];
            my ($exon_number) = grep { /^exon_number=/ } @featureInfo;
            my ($transcript_id) = grep { /^transcript_id=/ } @featureInfo;
            my ($transcript_type) = grep { /^transcript_type=/ } @featureInfo;

            if ($gene_name) {
                $gene_name =~ s/gene_name=//;
            }
            if ($exon_number) {
                $exon_number =~ s/exon_number=//;
            }
            if ($transcript_id) {
                $transcript_id =~ s/transcript_id=//;
            }
            if( $transcript_type) {
                $transcript_type =~ s/transcript_type=//;
            }

            my $outAnn = "gene_name=$gene_name;strand=$strand;transcript_id=$transcript_id;transcript_type=$transcript_type;exon_number=$exon_number";
            push @gffArr, $outAnn;
        }
        $gffAnnotation = join(";", @gffArr);
    }

    $tmp[7] = $tmp[7] . ";HG=$haplogroup;HS=$hsAnnotation;HP=$hpAnnotation;HV=$hvAnnotation;$gffAnnotation";
    $line = join("\t", @tmp);
    print OUT "$line\n";
}

close IN;
close OUT;

#################################

sub vcf2dict {
    my $input_vcf = shift;

    my %vcf_dict;

    open IN, " zcat $input_vcf | " or die "Cannot open $input_vcf: $!";
    while (my $line =<IN>) {
        chomp $line;
        next if $line =~ /^#/;  # Skip header line

        my @tmp = split("\t", $line);
        my $info = $tmp[7];

        my $coordinate = "chrM\t$tmp[1]\t$tmp[3]\t$tmp[4]";
        
        $vcf_dict{$coordinate} = $info;
    }
    close IN;

    return \%vcf_dict;
}


