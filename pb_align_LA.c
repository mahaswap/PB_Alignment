/*-------------------------------------------------------------

######### CODE WRITTEN FOR PB LOCAL ALIGNMENT BY MANOJ TYAGI and UPDATED BY SWAPNIL MAHAJAN #######

-------------------------------------------------------------*/


/*************************************************
	Program: pb_align_LA.c

Program for local PB sequence alignment, using the Smith-Waterman
algorithm and assuming a LINEAR gap penalty.
A traceback is used to determine the alignment, and
is determined by choosing that direction, for
which S(i-1,j-1)+sigma(a_i,b_j), S(i-1,j)+Delta and 
S(i,j-1)+Delta is maximum, i.e.  for which 

                    _
                   |
                   | H(i-1,j-1) + sigma(a_i,b_j)  (diagonal)
H(i,j) =  MAX of   | H(i-1,j)   + delta           (up)
                   | H(i,j-1)   + delta           (left)
                   | 0
                   | _


is a maximum.

*************************************************/


#include <stdio.h>
#include <ctype.h>   	// character handling
#include <stdlib.h>     // def of RAND_MAX 
#define N 1000         // size of protein arrays
#define AA 16           // number of amino acids
#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

// Agnel's PB Substitution Matrix
// Added by Swapnil Mahajan on 23-03-2012
const double sim[16][16] =
{{ 2.7832,-0.3529, 0.2985,-0.5713,-1.7361,-1.1270,-0.0934,-2.1397,-0.3089,-0.7020,-3.1815,-1.4108,-3.2729,-2.2225,-0.8006,-0.7240},
{ -0.3529, 2.7798,-0.7806,-1.2016,-0.8130,-2.2599,-0.6582, 0.0908, 0.7060,-0.7630,-0.4611,-0.0866,-3.8886,-1.5118,-1.0130,-0.2555},
{  0.2985,-0.7806, 2.0308,-0.3818,-1.6421,-1.1759, 0.5041,-2.0954,-1.4821,-1.4308,-3.1951,-3.8303,-2.7775,-1.8294,-0.6445,-0.0043},
{ -0.5713,-1.2016,-0.3818, 1.4084,-0.1811,-1.0630,-1.8802,-1.4975,-1.2259,-0.8866,-3.3702,-3.2946,-7.9873,-5.1749,-3.5603,-2.8899},
{ -1.7361,-0.8130,-1.6421,-0.1811, 3.0360, 0.7033, 0.8183, 0.6010,-2.0108,-0.2292,-0.7650,-2.8994,-5.8282,-2.2392,-3.0927,-3.5874},
{ -1.1270,-2.2599,-1.1759,-1.0630, 0.7033, 2.4402,-0.6799,-0.4917,-2.8870, 0.0011,-1.2239,-0.7720,-2.6212,-3.4147,-2.4299,-2.4166},
{ -0.0934,-0.6582, 0.5041,-1.8802, 0.8183,-0.6799, 3.2724,-1.0287,-0.4209,-1.1577,-2.3416,-0.9476,-0.8517, 0.3044, 0.4855, 1.0360},
{ -2.1397, 0.0908,-2.0954,-1.4975, 0.6010,-0.4917,-1.0287, 3.3588,-0.4706, 1.1012, 0.7138,-0.7945,-3.6712,-1.9999,-0.1331,-2.4846},
{ -0.3089, 0.7060,-1.4821,-1.2259,-2.0108,-2.8870,-0.4209,-0.4706, 3.6964, 1.0101,-0.6718,-0.0099,-3.9274,-2.0933,-1.5492, 0.7063},
{ -0.7020,-0.7630,-1.4308,-0.8866,-0.2292, 0.0011,-1.1577, 1.1012, 1.0101, 3.9565, 0.3985, 0.0186,-2.1680,-1.7781,-1.0124,-0.3773},
{ -3.1815,-0.4611,-3.1951,-3.3702,-0.7650,-1.2239,-2.3416, 0.7138,-0.6718, 0.3985, 2.7688,-0.6321,-1.6466,-1.6306,-3.0056,-2.0326},
{ -1.4108,-0.0866,-3.8303,-3.2946,-2.8994,-0.7720,-0.9476,-0.7945,-0.0099, 0.0186,-0.6321, 2.6433,-1.0826,-0.3362,-0.4381,-2.0649},
{ -3.2729,-3.8886,-2.7775,-7.9873,-5.8282,-2.6212,-0.8517,-3.6712,-3.9274,-2.1680,-1.6466,-1.0826, 1.1066,-0.5814,-1.3069,-1.1088},
{ -2.2225,-1.5118,-1.8294,-5.1749,-2.2392,-3.4147, 0.3044,-1.9999,-2.0933,-1.7781,-1.6306,-0.3362,-0.5814, 3.5846, 0.0375, 0.5285},
{ -0.8006,-1.0130,-0.6445,-3.5603,-3.0927,-2.4299, 0.4855,-0.1331,-1.5492,-1.0124,-3.0056,-0.4381,-1.3069, 0.0375, 3.4316,-0.0103},
{ -0.7240,-0.2555,-0.0043,-2.8899,-3.5874,-2.4166, 1.0360,-2.4846, 0.7063,-0.3773,-2.0326,-2.0649,-1.1088, 0.5285,-0.0103, 2.9327}};

/*
// Manoj's PB substitution Matrix
// matrix from december , bug removed
const double sim[16][16] =
{{2.28,-0.12,0.54,-0.29,-1.59,-0.54,0.31,-1.14,0.39,-1.15,-1.75,-0.60,-2.40,-1.40,-0.54,-0.36},
{-0.12,2.49,-0.21,-0.44,-0.48,-1.53,-0.73,0.20,0.24,0.32,-0.03,0.04,-2.98,-0.83,-0.55,0.33},
{0.54,-0.21,1.69,0.17,-1.10,-0.39,0.18,-1.63,-1.11,-1.03,-2.45,-2.21,-2.70,-1.68,-0.65,-0.01},
{-0.29,-0.44,0.17,1.35,-0.36,-0.49,-1.29,-1.20,-1.12,-0.92,-2.63,-1.56,-5.20,-3.07,-2.66,-2.10},
{-1.59,-0.48,-1.10,-0.36,3.05,0.75,1.37,0.66,-1.15,-0.76,-0.38,-1.76,-4.75,-0.58,-2.48,-2.22},
{-0.54,-1.53,-0.39,-0.49,0.75,2.21,-0.33,-0.34,-1.07,-0.34,-0.04,-0.33,-2.14,-1.99,-1.41,-1.91},
{0.31,-0.73,0.18,-1.29,1.37,-0.33,3.25,-0.74,-0.19,-0.51,-1.39,-0.74,-1.10,1.07,-0.01,0.47},
{-1.14,0.20,-1.63,-1.20,0.66,-0.34,-0.74,3.07,-0.92,1.18,0.51,-0.36,-2.93,-1.07,0.96,-1.81},
{0.39,0.24,-1.11,-1.12,-1.15,-1.07,-0.19,-0.92,3.37,1.54,-0.15,-0.22,-3.15,-0.97,-0.89,1.32},
{-1.15,0.32,-1.03,-0.92,-0.76,-0.34,-0.51,1.18,1.54,3.74,0.07,-0.12,-2.00,-0.44,-0.48,0.60},
{-1.75,-0.03,-2.45,-2.63,-0.38,-0.04,-1.39,0.51,-0.15,0.07,2.52,0.19,-1.02,-0.56,-1.71,-1.35},
{-0.60,0.04,-2.21,-1.56,-1.76,-0.33,-0.74,-0.36,-0.22,-0.12,0.19,2.24,-0.68,-0.27,0.06,-1.23},
{-2.40,-2.98,-2.70,-5.20,-4.75,-2.14,-1.10,-2.93,-3.15,-2.00,-1.02,-0.68,1.06,-0.77,-1.26,-1.10},
{-1.40,-0.83,-1.68,-3.07,-0.58,-1.99,1.07,-1.07,-0.97,-0.44,-0.56,-0.27,-0.77,3.65,0.26,0.36},
{-0.54,-0.55,-0.65,-2.66,-2.48,-1.41,-0.01,0.96,-0.89,-0.48,-1.71,0.06,-1.26,0.26,3.36,0.24},
{-0.36,0.33,-0.01,-2.10,-2.22,-1.91,0.47,-1.81,1.32,0.60,-1.35,-1.23,-1.10,0.36,0.24,2.83} };
*/

/*matrix with diagonals having double frequency */ 
/*const double sim[16][16] =
{{2.73,-0.38,0.26,-0.70,-1.87,-0.85,0.16,-1.43,0.14,-1.34,-2.08,-0.89,-2.91,-1.73,-0.85,-0.63},
{-0.38,2.91,-0.50,-0.87,-0.76,-1.86,-0.89,-0.10,-0.03,0.12,-0.36,-0.26,-3.50,-1.17,-0.88,0.04},
{0.26,-0.50,2.08,-0.27,-1.40,-0.73,0.01,-1.95,-1.39,-1.25,-2.81,-2.53,-3.24,-2.03,-0.99,-0.31},
{-0.70,-0.87,-0.27,1.47,-0.80,-0.97,-1.61,-1.65,-1.53,-1.28,-3.12,-2.01,-5.87,-3.56,-3.14,-2.54},
{-1.87,-0.76,-1.40,-0.80,3.45,0.42,1.19,0.34,-1.43,-0.97,-0.73,-2.08,-5.29,-0.93,-2.82,-2.52},
{-0.85,-1.86,-0.73,-0.97,0.42,2.53,-0.54,-0.69,-1.38,-0.59,-0.43,-0.69,-2.71,-2.37,-1.78,-2.25},
{0.16,-0.89,0.01,-1.61,1.19,-0.54,3.89,-0.93,-0.34,-0.60,-1.61,-0.93,-1.51,0.84,-0.23,0.30},
{-1.43,-0.10,-1.95,-1.65,0.34,-0.69,-0.93,3.42,-1.22,0.95,0.14,-0.70,-3.48,-1.43,0.61,-2.12},
{0.14,-0.03,-1.39,-1.53,-1.43,-1.38,-0.34,-1.22,3.81,1.35,-0.47,-0.51,-3.66,-1.30,-1.21,1.05},
{-1.34,0.12,-1.25,-1.28,-0.97,-0.59,-0.60,0.95,1.35,4.30,-0.19,-0.36,-2.45,-0.71,-0.74,0.39},
{-2.08,-0.36,-2.81,-3.12,-0.73,-0.43,-1.61,0.14,-0.47,-0.19,2.81,-0.19,-1.61,-0.96,-2.10,-1.70},
{-0.89,-0.26,-2.53,-2.01,-2.08,-0.69,-0.93,-0.70,-0.51,-0.36,-0.19,2.60,-1.24,-0.64,-0.29,-1.55},
{-2.91,-3.50,-3.24,-5.87,-5.29,-2.71,-1.51,-3.48,-3.66,-2.45,-1.61,-1.24,0.98,-1.35,-1.83,-1.64},
{-1.73,-1.17,-2.03,-3.56,-0.93,-2.37,0.84,-1.43,-1.30,-0.71,-0.96,-0.64,-1.35,3.94,-0.13,0.01},
{-0.85,-0.88,-0.99,-3.14,-2.82,-1.78,-0.23,0.61,-1.21,-0.74,-2.10,-0.29,-1.83,-0.13,3.67,-0.10},
{-0.63,0.04,-0.31,-2.54,-2.52,-2.25,0.30,-2.12,1.05,0.39,-1.70,-1.55,-1.64,0.01,-0.10,3.23} }; */

main(int argc, char *argv[]) {
	// function prototypes
	void error(char *);		/** error handling */
	int char2AA(char);
	char AA2char(int);

	// variable declarations
	FILE * in1, *in2, *pam;
	//FILE * output;
	char ch;
	int temp;
	int i,j,tempi,tempj,x,y;
	float diag,down,right,DELTA;
	int topskip, bottomskip;
	char aout[N],bout[N];
	int Aend,Bend,Abegin,Bbegin;
	int xMax, yMax;
	double max, Max;
	int pblen1, pblen2;
		// Max is first found maximum in similarity matrix H
		// max is auxilliary to determine largest of
		// diag,down,right, xMax,yMax are h-coord of Max
	//short a[N], b[N];
	char *a1, *b1;
	char a[N], b[N];
	float h[N][N];
	//int sim[AA][AA];		// PAM similarity matrix
	short xTraceback[N][N], yTraceback[N][N];

	/**** Error handling for input file ****/
	if (argc != 6) {
	     	fprintf(stderr,"%s PBseq1 PBseq2 len1 len2 gapPenalty\n",argv[0]);
		exit(1);
		}

	/***** Initialization of input file and pair array **/
	//output = fopen(argv[6], "w");
	/* assign PB sequences to "a" & "b" */
	a1 = argv[1];
	b1 = argv[2];
	pblen1 = atoi(argv[3]);
	pblen2 = atoi(argv[4]);
	DELTA = atof(argv[5]);

	
	/*printf("i am here\n");*/
	for(i=1;i<=pblen1;i++) {
		a[i]=a1[i-1];
		/*printf("%c",a[i]);*/
		}
	for(i=1;i<=pblen2;i++) {
		b[i]=b1[i-1];
		}

	/** initialize traceback array **/
	Max=xMax=yMax=0;
	for (i=0;i<=pblen1;i++)
		for (j=0;j<=pblen2;j++) {
			xTraceback[i][j]=-1;
			yTraceback[i][j]=-1;
			}

	/** compute "h" local similarity array **/
	for (i=0;i<=pblen1;i++) h[i][0]=0;
	for (j=0;j<=pblen2;j++) h[0][j]=0;
	
	for (i=1;i<=pblen1;i++)
		for (j=1;j<=pblen2;j++) {
			diag    = h[i-1][j-1] + sim[a[i]-65][b[j]-65];
			down    = h[i-1][j] + DELTA;
			right   = h[i][j-1] + DELTA;
			max=MAX3(diag,down,right);
			if (max <= 0)  {
				h[i][j]=0;
				xTraceback[i][j]=-1;
				yTraceback[i][j]=-1;
					// these values already -1
				}
			else if (max == diag) {
				h[i][j]=diag;
				xTraceback[i][j]=i-1;
				yTraceback[i][j]=j-1;
				}
			else if (max == down) {
				h[i][j]=down;
				xTraceback[i][j]=i-1;
				yTraceback[i][j]=j;
				}
			else  {
				h[i][j]=right;
				xTraceback[i][j]=i;
				yTraceback[i][j]=j-1;
				}
			if (max > Max){
				Max=max;
				xMax=i;
				yMax=j;
				}
			}  // end for loop


	// output values for gnuplot
	i=xMax; j=yMax;
	/*while ( i>0 && j>0 && h[i][j]>0){
		printf("%d %d\n",i,j);
		//printf("%c %c\n",AA2char(a[i]),AA2char(b[j]));
		printf("score %f\n",h[i][j]);
			// previous 2 lines for debugging

		tempi=i;
		tempj=j;
		i=xTraceback[tempi][tempj];
		j=yTraceback[tempi][tempj];
			// WARNING -- following 2 lines incorrect!
			// You need tempi, tempj
			//i=xTraceback[i][j];
			//j=yTraceback[i][j];
	}*/




	// initialize output arrays to be empty -- this is unnecessary
	/*for (i=0;i<N;i++) aout[i]=' ';
	for (i=0;i<N;i++) bout[i]=' ';*/


	// reset to max point to do alignment
	i=xMax; j=yMax;
	x=y=0;
	topskip = bottomskip = 1;
	while (i>0 && j>0 && h[i][j] > 0){
		topskip    = (j>yTraceback[i][j]);
		bottomskip = (i>xTraceback[i][j]);
		if (topskip && bottomskip) {
			/*aout[x++]=AA2char(a[i]);
			bout[y++]=AA2char(b[j]);*/
			aout[x++]=(a[i]);
			bout[y++]=(b[j]);
			}
		else if (topskip) {
			aout[x++]='-';
			/*bout[y++]=AA2char(b[j]);*/
			bout[y++]=(b[j]);
			}
		else if (bottomskip) {
			/*aout[x++]=AA2char(a[i]);*/
			aout[x++]=(a[i]);
			bout[y++]='-';
			}
		
		tempi=i;tempj=j;
		i=xTraceback[tempi][tempj];
		j=yTraceback[tempi][tempj];
			// Warning -- following 2 lines no good
			// i=xTraceback[i][j];
			// j=yTraceback[i][j];
	}



	// print alignment
	/*printf("\n");
	printf("\n");*/
	printf("\n");
	// print alignment
	for (i=x-1;i>=0;i--) 
	{
	  //fprintf(output,"%c",aout[i]); 
	  printf("%c",aout[i]);
	}
	//fprintf(output,"\n");
	printf("\n");
	for (j=y-1;j>=0;j--) 
	{
	  //fprintf(output,"%c",bout[j]); 
	  printf("%c",bout[j]);
	}
	//fprintf(output,"\n");
	printf("\n");
	/*printf("\n");*/
	printf("raw score:%2.3f:%2.3f:%d\n",Max,(Max/x),x);

}

void error(char * s) {
	fprintf(stderr,"%s\n",s);
	exit(1);
}



int char2AA(char ch){
        switch (ch) {
                case 'A' : return 0;
                case 'B' : return 1;
                case 'C' : return 2;
                case 'D' : return 3;
                case 'E' : return 4;
                case 'F' : return 5;
                case 'G' : return 6;
                case 'H' : return 7;
                case 'I' : return 8;
                case 'J' : return 9;
                case 'K' : return 10;
                case 'L' : return 11;
                case 'M' : return 12;
                case 'N' : return 13;
                case 'O' : return 14;
                case 'P' : return 15;
                default  :
                        fprintf(stderr,"Nonstandard PROTEIN BLOCK code: %d\n",ch);
                        exit(1);
                }
}


char AA2char(int x){
        switch (x) {
                case 0 : return 'A';
                case 1 : return 'B';
                case 2 : return 'C';
                case 3 : return 'D';
                case 4 : return 'E';
                case 5 : return 'F';
                case 6 : return 'G';
                case 7 : return 'H';
                case 8 : return 'I';
                case 9 : return 'J';
                case 10: return 'K';
                case 11: return 'L';
                case 12: return 'M';
                case 13: return 'N';
                case 14: return 'O';
                case 15: return 'P';
                default  : fprintf(stderr,"Bad char: %d\n",x);
                           error("Bad integer representation of AA");
                }
}




/*---------------------------------------------------
Output of
	a.out Cdc25 Ste5 pam250 -10


999 821
998 820
997 819
996 818
995 817
994 816
993 815
992 814
991 813
990 812
989 811
988 810
987 809
986 808
985 807
984 806
983 805
982 804
981 803
980 802
979 801
978 800
977 799
976 798
975 797
974 796
973 795
972 794
971 793
970 792
969 791
968 790
967 789
966 788
965 787
964 786
963 785
962 784
961 783
960 782
959 781
958 780
957 779
956 778
955 777
954 776
953 775
952 774
951 773
950 772
949 771
948 770
947 769
946 768
945 767
944 766
943 765
943 764
942 763
941 762
941 761
940 760
939 759
938 758
937 757
936 756
935 755
934 754
933 753
932 752
931 751
930 750
929 749
928 748
927 747
926 746
925 745
924 744
923 743
922 742
921 741
920 740
919 739
918 738
917 737
916 736
915 735
914 734
913 733
912 732
911 731
910 730
909 729
908 728
908 727
907 726
906 725
905 724
904 723
903 722
902 721
901 720
900 719
899 718
898 717
897 716
896 715
895 714
894 713
893 712
892 711
891 710
890 709
889 708
888 707
887 706
886 705
885 704
884 703
883 702
882 701
881 700
880 699
879 698
878 697
877 696
876 695
875 694
874 693
873 692
872 691
871 690
870 689
869 688
868 687
867 686
866 685
865 684
864 683
863 682
862 681
861 680
860 679
859 678
858 677
857 676
856 675
855 674
854 673
853 672
853 671
852 670
851 669
850 668
849 667
848 666
847 665
846 664
845 663
844 662
843 661
842 660
841 659
840 658
839 657
838 656


EHLKIISKPKSRIRN-LEINSSTYEQINQNVLLLEILENLDLSIFINLKNLIKTPSILLDLESEEFLVHAM-SSVSSVLTEFFDIKQAFHDIVIRLIMTTQQTTL-DD-PYLFSSMRSNFPVGHHEPFKNISNTPLVKGPFHKKNEQLALSLFHVLVSQDVEFNNL
DELVLLLPPREIAKQLCILEFQSFSHISRIQFLTKIWDNLNRFSPKEKTSTFYLSNHLVNFVTETIVQEEEPRRRTNVLAYFIQVCDYLRELNNFASLFSIISALNSSPIHRLRKTWANLNSKTLASFELLNNLTEARKNFSNYRDCLENCVLPCVPFLGVYFTDL

---------------------------------------------------*/
