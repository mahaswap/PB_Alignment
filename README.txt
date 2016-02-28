##################################################
Written by Swapnil MAHAJAN. 26th March 2013.
##################################################

## Compile C scripts for PB-ALIGN (You can use "cc" or "gcc" compiler on linux and Mac-OS; for Windows OS you will need visual basic suite to compile these C scripts)

gcc -g -o pb_align_GA pb_align_GA.c -lm
gcc -g -o pb_align_LA pb_align_LA.c -lm
##################################################

##PB-ALIGN pairwise global alignment

perl pb_align_pairwise_GA.pl path/PB_seq_file1 path/PB_seq_file2 path/output_file_name
e.g.
perl pb_align_pairwise_GA.pl ./d1aq0a_.pbseq ./d1aq1a_.pbseq ./d1aq0a_.aln
##################################################

##PB-ALIGN pairwise local alignment

perl pb_align_pairwise_LA.pl path/PB_seq_file1 path/PB_seq_file2 path/output_file_name
e.g.
perl pb_align_pairwise_LA.pl ./d1aq0a_.pbseq ./d1aq1a_.pbseq ./d1aq0a_.aln
##################################################

##PB-ALIGN database mining using global alignment
# You need to install Statistics::Basic Perl module to use these following scripts
# You can have one or more input files in input folder. But all the input files should be PB sequence files with .pbseq extesion
# astral95_1.75A folder (or a zipped file. Please unzip it before using) is provided which contains PB seq and AA seq files for Astral95 dataset from SCOP 1.75A.
# Last two arguments for length_cutoff and ranks are optional, default values are set to 30 and 20 respectively.

perl pb_align_mining_GA.pl path/input_PB_seq_folder ./astral95_1.75A/PB_SEQ path/output_file_name length_cutoff[optional, default 30] ranks[optional, default 20]

#Similarly for database mining using locl alignment
perl pb_align_mining_LA.pl path/input_PB_seq_folder ./astral95_1.75A/PB_SEQ path/output_file_name length_cutoff[optional, default 30] ranks[optional, default 20]

#Output file format for these two scripts are same.
#Two output files are created output_file.sc and output_file.aln.
#These two output files have fields seperated by "#" symbol.

#**** Normalized alignment score cutoff is -0.25 for significant alignments.(Tyagi et al, Proteins 2007) ****#

### output_file.sc format

input_file_id
input_seq_length
raw_alignment_score
normalized_alignment_score
first_hits_family_id
Z_score
total_number_of_alignments_performed


### output_file.aln format

input_file_id
input_seq_length
hits_file_id
hits_seq_length
hits_family_id
raw_alignment_score
normalized_alignment_score
alignment_length
aligned_input_seq
aligned_hits_seq
Z_score
total_number_of_alignments_performed

