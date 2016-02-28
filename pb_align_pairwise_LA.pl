#!/usr/bin/perl 
=head1
Written by Swapnil MAHAJAN. 26th March 2013.
=cut

use File::Basename;
use List::Util qw(first);
#use Statistics::Basic::Mean;
#use Statistics::Basic::StdDev;
use threads;
use threads::shared;

my ($aa_file):shared=$ARGV[0];  # input directory with .pbseq files
my ($pb_file):shared=$ARGV[1];  # ast95_1.75A directory for PB sequences
my ($out_file):shared=$ARGV[2]; # output file name
if($out_file eq "")
{$out_file='./pb_align_pairwise_LA_result';}

my ($GAP):shared= -5.0;
#my ($id_file):shared='./astral95_1.75A/scop_1.75A_domain_new.id'; # file with fold ids
#my ($RANKS):shared=20;  # top ranks to be displayed
#my $ext='.pbseq';  # extension of input files
#my $THREADS=4; # number of threads to be run simultaneously

open(OUT_SC, ">$out_file")||die("ERR:$0:OUT:$!\n");
#open(OUT_ALN, ">>$out_file\.aln")||die("ERR:$0:OUT:$!\n");
#opendir(PB_DIR, $pb_dir)||die("ERR:$0:PB_DIR:$!\n");

my $aa_id="";my $fold="";

( $aa_id, $path, $suffix) = fileparse($aa_file, qr/\.[^.]*$/);
undef $path; undef $suffix;

#$fold=`grep -m 1 $aa_id $id_file| cut -c 9-`;
#chomp $fold;
#$fold=~s/ //g;
#$fold=~s/\d+\.\d+$//;

$aa_seq=""; $aa_len=0;
open(AA, "$aa_file")||die("ERR:$0:AA:$!\n");
foreach my $l(<AA>)
{
  if($l!~/^>/){chomp $l; $aa_seq.=$l;}
}
close (AA);
$aa_seq=~s/[Z|X]//gi;
$aa_len=length($aa_seq);

$pb_len=0;

  $pb_seq="";
  open(PB, "$pb_file")||die("ERR:$0:PB:$!\n");
  foreach my $l(<PB>)
  {
    if($l!~/^>/){chomp $l; $pb_seq.=$l;}
  }
  close(PB);
  $pb_seq=~s/[Z|X]//gi;
  $pb_len=length($pb_seq);
  
   @aln=();
  #@aln=`./threading $aa_file $pb_dir/$pb_file $GAP`;
  @aln=`./pb_align_LA $aa_seq $pb_seq $aa_len $pb_len $GAP`;
  chomp @aln;
  print OUT_SC "$aln[0]\n";
  print OUT_SC "$aln[1]\n";
  print OUT_SC "$aln[2]\n";
  print OUT_SC "$aln[3]\n";

