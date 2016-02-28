#!/usr/bin/perl 
=head1
Written by Swapnil MAHAJAN. 26th March 2013.
=cut

use File::Basename;
use List::Util qw(first);
use Statistics::Basic::Mean;
use Statistics::Basic::StdDev;
use threads;
use threads::shared;

my ($aa_dir):shared=$ARGV[0];  # input directory with .pbseq files
my ($pb_dir):shared=$ARGV[1];  # ast95_1.75A directory for PB sequences
my ($out_file):shared=$ARGV[2]; # output file name
if($out_file eq "")
{$out_file='./pb_align_GA_result';}

my ($LEN_THRESHOLD):shared=$ARGV[3];
chomp $LEN_THRESHOLD;
if ($LEN_THRESHOLD eq "")
{$LEN_THRESHOLD=30;}

my ($RANKS):shared= $ARGV[4];  # top ranks to be displayed
if ($RANKS eq "")
{$RANKS = 20;}

my ($GAP):shared= -3.0;
my ($id_file):shared='./astral95_1.75A/scop_1.75A_domain_new.id'; # file with fold ids
my $ext='.pbseq';  # extension of input files
my $THREADS=4; # number of threads to be run simultaneously

open(OUT_SC, ">>$out_file\.sc")||die("ERR:$0:OUT:$!\n");
open(OUT_ALN, ">>$out_file\.aln")||die("ERR:$0:OUT:$!\n");
#opendir(PB_DIR, $pb_dir)||die("ERR:$0:PB_DIR:$!\n");
opendir(AA_DIR, $aa_dir)||die("ERR:$0:AA_DIR:$!\n");

my @aaseq_files=(grep /$ext/,readdir(AA_DIR));
closedir(AA_DIR);

for($i=0;$i<@aaseq_files;$i=$i+$THREADS)
{
  my @threads=();
  #push(@threads, threads->new(\&pb_align,$aaseq_files[$i]));
  for($j=0;$j<$THREADS;$j++)
  {
    if(($i+$j) <@aaseq_files){push(@threads, threads->new(\&pb_align,$aaseq_files[$i+$j]));}
  }
  
  foreach my $myThread(@threads)
  {
    ($out_sc, $out_aln) = $myThread->join();
    sleep(1);
    print OUT_SC "$out_sc\n";
    foreach my $l(@{$out_aln}){print OUT_ALN "$l\n";}
  }
}


sub pb_align{

my $aa_file=$_[0];
my $aa_id="";my $fold="";

( $aa_id, $path, $suffix) = fileparse($aa_file, qr/\.[^.]*$/);
undef $path; undef $suffix;

#$fold=`grep -m 1 $aa_id $id_file| cut -c 9-`;
#chomp $fold;
#$fold=~s/ //g;
#$fold=~s/\d+\.\d+$//;

$aa_seq=""; $aa_len=0;
open(AA, "$aa_dir/$aa_file")||die("ERR:$0:AA:$!\n");
foreach my $l(<AA>)
{
  if($l!~/^>/){chomp $l; $aa_seq.=$l;}
}
close (AA);
$aa_seq=~s/[Z|X]//gi;
$aa_len=length($aa_seq);

@score=(); @id=(); @sc_idx=(); @aln_sc=(); @norm_sc=();
$cnt_aln=0;$pb_len=0;

opendir(PB_DIR, $pb_dir)||die("ERR:$0:PB_DIR:$!\n");
foreach $pb_file(grep /\.pbseq$/,readdir(PB_DIR))
{
  $pb_seq="";
  open(PB, "$pb_dir/$pb_file")||die("ERR:$0:PB:$!\n");
  foreach my $l(<PB>)
  {
    if($l!~/^>/){chomp $l; $pb_seq.=$l;}
  }
  close(PB);
  $pb_seq=~s/[Z|X]//gi;
  $pb_len=length($pb_seq);
  
  $percent=0;
  $percent=(100-((abs($aa_len-$pb_len)/$aa_len)*100));
  if($percent < $LEN_THRESHOLD){next;}
  
  @aln=();
  #@aln=`./threading $aa_file $pb_dir/$pb_file $GAP`;
  @aln=`./pb_align_GA $aa_seq $pb_seq $aa_len $pb_len $GAP`;
  if(@aln)
  {
    $cnt_aln++;
    chomp @aln;
    $sc=0;$n_sc=0;$aln_len=0;$pb_id="";
    @tmp=();
    @tmp=split(/\:/,$aln[3]);
    $sc=$tmp[1]; $n_sc=$tmp[2]; $aln_len=$tmp[3];
    $sc= sprintf("%.3f",$sc); $n_sc= sprintf("%.3f",$n_sc);
    #$z_val= $tmp[4];

    push(@score,$sc);
    push(@norm_sc,$n_sc);
    #push(@zval,$z_val);
    ( $pb_id, $path, $suffix) = fileparse($pb_file, qr/\.[^.]*$/);
    undef $path; undef $suffix;
    $pb_fold="";
    $pb_fold=`grep -m 1 $pb_id $id_file| cut -c 9-`;
    chomp $pb_fold;
    $pb_fold=~s/ //g;
    #$pb_fold=~s/\d+\.\d+$//;
    push(@id,$pb_fold);

    push(@aln_sc,"$aa_id\#$aa_len\#$pb_id\#$pb_len\#$pb_fold\#$sc\#$n_sc\#$aln_len\#$aln[1]\#$aln[2]");
  }
}
closedir(PB_DIR);

@sc_idx =  sort { $score[$b] <=> $score[$a] } 0..$#score;
#@zval  = @zval [ @sc_idx ];
@score  = @score [ @sc_idx ];
@norm_sc  = @norm_sc [ @sc_idx ];
@id  = @id [ @sc_idx ];
@aln_sc  = @aln_sc [ @sc_idx ];
#$hit_rank = first { $id[$_] eq $fold } 0..$#id;
#$hit_rank+=1;

my $mean=new Statistics::Basic::Mean(\@score)->query;
my $std_dev=new Statistics::Basic::StdDev(\@score)->query;
$zval=sprintf("%.2f",(($score[0]-$mean)/$std_dev));

my $out_sc="$aa_id\#$aa_len\#$score[0]\#$norm_sc[0]\#$id[0]\#$zval\#$cnt_aln";
my @out_aln=();
for (my $i=0;$i<$RANKS;$i++)
{
  $zval="";
  @tp=();
  @tp=split(/\#/,$aln_sc[$i]);
  $zval=sprintf("%.2f",(($tp[5]-$mean)/$std_dev));
  
  $out_aln[$i]="$aln_sc[$i]"."\#$zval\#$cnt_aln";
}
return($out_sc,\@out_aln);
}
