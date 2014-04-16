#after format transfered, the tag, like barcode, at 5` should be cut-off. 
if(@ARGV!=2){
	die "usage : plesse input the file with format 'FASTQ' and input the cut-off length\n
example : perl filter_fastq_bam.pl file1 30; #A read with length < 30 will be filtered.\n
author: Yang Li\n
please feel free to contact me: yeli7068\@gmail.com\n";
	exit;}
open IN,"$ARGV[0]";
open OUT,">$ARGV[0].up_qc";
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
					
		s/\A(.){$ARGV[1]}//;
		$len_seq = length($_);
		print OUT "$_\n";
	}elsif ($n == 3) {
		print OUT "$_\n";
		$n = 4;
	}elsif ($n == 4) {
		$n = 1;
		$len_qua = length ($_);
		if ($len_seq == $len_qua) {
			print OUT "$_\n";
		} else {
			$len_cutoff = $len_qua - $len_seq;
			#print "cufoff is $len_cutoff\nlen_qua is $len_qua\nlen_seq is $len_seq\n";
			s/\A(.){$len_cutoff}//;  #将质量值补齐
			print OUT "$_\n";
		}
	}
}
