#after tags being removed, some of reads might be empty or shorter length for downstream analysis.
#this program will help you to cut off the reads with length lower than you wished.
if(@ARGV!=2){
	die "usage : plesse input the file with format 'FASTQ' and input the cut-off length\n
example : perl filter_fastq_bam.pl file1 30; #A read with length < 30 will be filtered.\n
author: Yang Li\n
please feel free to contact me: yeli7068\@gmail.com\n";
	exit;}
open IN,"$ARGV[0]";
open OUT,">$ARGV[0].length$ARGV[1]_filter_fq";
$fq = "";
$n = 1;
while(<IN>){
  chomp;
	if ($n == 1) {
		$n = 2;
		$fq .= "$_\n";}
	elsif ($n == 2) {
		$len = length($_);
		$n = 3;
		$fq .= "$_\n";}
	elsif ($n == 3) {
		$n = 4;
		$fq .= "$_\n";}
	elsif ($n == 4) {
		$fq .= "$_\n";
		$n = 1;
		if($len >= $ARGV[1]){
			print OUT "$fq";
			$fq = "";
			
		}else {
			$fq = "";
			
		}
	} 
}

close IN;
close OUT;
