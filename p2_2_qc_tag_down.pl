#this program was developed to cutoff the nucl at 3`. if set 270, part of the read which was longer than 270 will be cut-off. 
if(@ARGV!=2){
	die "usage : plesse input the file with format 'FASTQ' and input the cut-off length\n
example : perl filter_fastq_bam.pl file1 30; #A read with length < 30 will be filtered.\n
author: Yang Li\n
please feel free to contact me: yeli7068\@gmail.com\n";
	exit;}
open IN,"$ARGV[0]";
open OUT,">$ARGV[0].down_qc";
$n = 1;
while (<IN>) {
	chomp;
	if ($n == 1) {
		print OUT "$_\n";
		$n = 2;
	}elsif ($n == 2) {
		$n = 3;
		#s/\A(T)*(A)*G(C)*(G)*AGCTCTGCAGATATC(A|G|C|T){10}//;
		#s/\A(T)*(A)*//;
		#$seq = $_;
		$len_seq = length ($_);		
		if ($len_seq > $ARGV[1]) {
			$cutoff_down = $len_seq - $ARGV[1];
			s/(.){$cutoff_down}$//;
			}
		print OUT "$_\n";
		$len_after_cutoff = length ($_);
	}elsif ($n == 3) {
		print OUT "$_\n";
		$n = 4;
	}elsif ($n == 4) {
		$n = 1;
		$len_qua = length ($_);
		if ($len_after_cutoff == $len_qua) {
			print OUT "$_\n";
		} else {
			$len_cutoff = $len_qua - $len_after_cutoff;
			#print "cufoff is $len_cutoff\nlen_qua is $len_qua\nlen_seq is $len_seq\n";
			s/(.){$len_cutoff}$//;  #将质量值补齐
			print OUT "$_\n";
		}
	}
}
