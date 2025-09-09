#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

$ENV{LC_ALL} = 'C';

# ---------------- CLI options & defaults ----------------
my $ref        = '/media/gencardio/KINGSTON/hg38/hg38.fa';
my $contig     = 'chrM';                               # e.g. MT or chrM (omit if ref is chrM-only)
my $outdir     = '/media/gencardio/KINGSTON/Mitochondrial_simulation';
my $coverages  = '5000,1000,500,100,30,10,5,1,0.1';
my $readlen    = 100;                              # ART read length (R1/R2)
my $ins_mean   = 170;                              # ART -m (fragment mean)
my $ins_sd     = 60;                               # ART -s (fragment stdev)
my $system     = 'HS25';                           # ART system (HiSeq 2500): HS25, HS20, MSv3, NovaSeq...
my $n_snvs     = 25;                               # number of random SNVs to generate
my $seed_base  = 2025;     #Set 1 43  Set 2 2025   # base RNG seed for SNV generation
my $prefix     = 'RB';                               # output filename prefix (defaults to contig id)
my $help       = 0;

GetOptions(
  'ref=s'       => \$ref,
  'contig=s'    => \$contig,
  'outdir=s'    => \$outdir,
  'coverages=s' => \$coverages,
  'readlen=i'   => \$readlen,
  'ins=i'       => \$ins_mean,
  'sd=i'        => \$ins_sd,
  'n-snvs=i'    => \$n_snvs,
  'prefix=s'    => \$prefix,
  'help|h'      => \$help,
) or die "Bad options. Use --help\n";

if ($help || $ref eq '') {
  print <<"USAGE";
Usage:
  $0 --ref <reference.fasta> [--contig MT|chrM] [options]

Examples:
  $0 --ref GRCh38.fa --contig MT
  $0 --ref chrM.fa --n-snvs 40 --readlen 101 --system NovaSeq --outdir out

Options:
  --outdir      Output directory (default: $outdir)
  --coverages   Comma list (default: $coverages)
  --readlen     Read length for R1/R2 (default: $readlen)
  --ins         Fragment mean (default: $ins_mean)
  --sd          Fragment stdev (default: $ins_sd)
  --n-snvs      Number of random SNVs (default: $n_snvs)
  --prefix      Output filename prefix (default: contig id)
USAGE
  exit 0;
}

# ---------------- tiny helpers ----------------
sub which {
  my ($exe) = @_;
  my $rc = system("which $exe >/dev/null 2>&1");
  if ($rc != 0) { die "Missing dependency: $exe\n"; }
}

sub slurp_cmd {
  my ($cmd) = @_;
  my $out = `$cmd`;
  if ($? != 0) { die "Failed: $cmd\n"; }
  chomp $out;
  return $out;
}

# Generate a simple VCF of N random SNVs on chrM (avoid termini).
# Uses REF base from fasta; picks ALT != REF.

sub generate_random_snvs {
  my ($mt_fa, $mt_id, $mt_len, $n, $seed, $vcf_out) = @_;

  open(my $V, '>', $vcf_out) or die "Cannot write $vcf_out: $!\n";
  print $V "##fileformat=VCFv4.2\n";
  print $V "##contig=<ID=$mt_id,length=$mt_len>\n";
  print $V "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

  my %seen;
  my @bases = ('A','C','G','T');
  srand($seed);
  my $made = 0;

  while ($made < $n) {
    my $pos = 2 + int(rand($mt_len - 2)); # skip 1 and last
    if (exists $seen{$pos}) { next; }
    $seen{$pos} = 1;

    my $cmd = "samtools faidx '$mt_fa' '$mt_id:$pos-$pos' | tail -n +2 | tr -d '\\n' | tr 'acgt' 'ACGT'";
    my $refb = slurp_cmd($cmd);
    if ($refb !~ /^[ACGT]$/) { next; }

    my @alts = ();
    for my $b (@bases) {
      if ($b ne $refb) { push @alts, $b; }
    }
    my $altb = $alts[int(rand(@alts))];

    print $V "$mt_id\t$pos\t.\t$refb\t$altb\t.\t.\t.\n";
    $made++;
  }
  close $V;
  print " INFO: SNV VCF : $vcf_out (N=$n)\n";
}

# ---------------- check tools ----------------
which('samtools');
which('bcftools');
which('art_illumina');
which('gzip');
which('grep');
which('sed');

# make output dir
my $mk = "mkdir -p '$outdir'";
system($mk);

# ---------------- extract chrM fasta ----------------
my $mt_fa  = "$outdir/chrM.fa";
my $mt_id  = $contig;

if ($contig ne '') {
  my $cmd = "samtools faidx '$ref' '$contig' > '$mt_fa'";
  print " INFO: $cmd\n";
  my $rc = system($cmd);
  if ($rc != 0) { die "Could not extract contig $contig from $ref\n"; }
  $mt_id = $contig;
} 
else {
  my $n_hdr = slurp_cmd("grep -c '^>' '$ref'");
  if ($n_hdr > 1) {
    my $found = 0;
    for my $try ('MT','chrM','M') {
      my $cmd = "samtools faidx '$ref' '$try' > '$mt_fa'";
      print " INFO: $cmd\n";
      my $rc = system($cmd);
      if ($rc == 0) { 
        $mt_id = $try; 
        $found = 1; 
        last; 
      }
    }
    if ($found == 0) { 
      print " ERROR: Reference has multiple contigs; please pass --contig (e.g., MT or chrM)\n"; 
      exit;
    }
  } 
  else {
    my $cp = "cp '$ref' '$mt_fa'";
    print " INFO: $cp\n"; 
    system($cp) == 0 or die "Failed: $cp\n";
    $mt_id = slurp_cmd("grep -m1 '^>' '$mt_fa' | sed 's/^>//; s/[[:space:]].*\$//'");
  }
}

my $cmd = "samtools faidx '$mt_fa'";
system($cmd);

my $mt_len = slurp_cmd("cut -f2 '$mt_fa.fai'");

if ($prefix eq '') { 
    $prefix = $mt_id; 
}
print " INFO: chrM: $mt_id  len=$mt_len bp  prefix=$prefix  outdir=$outdir\n";

# ---------------- generate random SNVs (VCF) ----------------

my $vcf_raw = "$outdir/chrM.mutations.vcf";
generate_random_snvs($mt_fa, $mt_id, $mt_len, $n_snvs, $seed_base, $vcf_raw);

# compress + index VCF
my $vcf_gz = "$vcf_raw.gz";
$cmd  = "bcftools sort -Oz -o '$vcf_gz' '$vcf_raw' > /dev/null";
system($cmd);

$cmd = "bcftools index -t '$vcf_gz' > /dev/null";
system($cmd);

# ---------------- make mutated chrM fasta ----------------

print " INFO: Creating mutated fasta\n";
my $mut_fa = "$outdir/chrM.mut.fa";
$cmd   = "bcftools consensus -f '$mt_fa' '$vcf_gz' -o '$mut_fa' > /dev/null";
system($cmd);
print " INFO: $cmd\n";
print " INFO: Mutated fasta : $mut_fa\n";

$cmd = "samtools faidx '$mut_fa'";
print " INFO: $cmd\n"; 
system($cmd);

# ---------------- parse coverages ----------------

my @covs = ();
my @raw = split(/,/, $coverages);
for my $c (@raw) {
  $c =~ s/^\s+//;
  $c =~ s/\s+$//;
  if ($c ne '') {
    push @covs, $c + 0;
  }
}
if (scalar(@covs) == 0) { 
  print " ERROR: No valid coverages parsed from --coverages\n"; 
  exit;
}

# ---------------- simulate with ART ----------------

my $i = 0;
for my $cov (@covs) {
  if (!($cov > 0)) { 
    print " ERROR: Coverage must be positive: $cov\n";
    exit;
  }

  my $tag = $cov;
  $tag =~ s/\.0$//;

  my $art_prefix = "$outdir/${prefix}${tag}_s4_";

  # my $art_prefix = "$outdir/${prefix}.mut.art.${tag}X";
  print " INFO: ART: cov=$cov× len=$readlen m=$ins_mean s=$ins_sd system=$system\n";

  my $art_cmd = "art_illumina -ss $system -i '$mut_fa' -p -l $readlen -f $cov -m $ins_mean -s $ins_sd -o '$art_prefix' > /dev/null";
  print " INFO: $art_cmd\n"; 
  system($art_cmd);;

  # remove .aln generated by ART (if present)
  my $aln = "${art_prefix}.aln";
  my $rm_aln = "rm -f '$aln'";
  print " INFO: $rm_aln\n"; 
  system($rm_aln);

  # rename and gzip FASTQs
  my $r1_src = "${art_prefix}1.fq";
  my $r2_src = "${art_prefix}2.fq";
  
  my $r1     = "$outdir/${prefix}${tag}_s1_R1.fastq";
  my $r2     = "$outdir/${prefix}${tag}_s1_R2.fastq";


  my $mv1 = "mv '$r1_src' '$r1'";
  print " INFO: $mv1\n"; 
  system($mv1);

  my $mv2 = "mv '$r2_src' '$r2'";
  print " INFO: $mv2\n"; 
  system($mv2);

  my $gz = "gzip -f '$r1' '$r2'";
  print " INFO: $gz\n"; 
  system($gz);

  # create subfolder for this coverage and move gz files inside
  my $covdir = "$outdir/${tag}X";
  my $mkcov  = "mkdir -p '$covdir'";
  print " INFO: $mkcov\n"; 
  system($mkcov);

  my $mv_gz1 = "mv '$r1.gz' '$covdir/'";
  print " INFO: $mv_gz1\n"; 
  system($mv_gz1);

  my $mv_gz2 = "mv '$r2.gz' '$covdir/'";
  print " INFO: $mv_gz2\n";
  system($mv_gz2);

  $i = $i + 1;
}

# remove .aln files
my $rm_all_aln = "rm -f '$outdir'/*.aln";
print " INFO: $rm_all_aln\n"; 
system($rm_all_aln);
