#afer p2_1 2 3, the data at 3` end might still be unsatisfied. This program was to cut off the nucl at this 
#end with low quality (like 10).  This program will detect length as user defined from 3 end. And it will check each
#nucl with quality. If a nucl with quality lower than 10, the program will cutoff this nucl and the next nucls.
use strict;
#use warnings;
if (@ARGV != 2) {
	die 
"usage: perl qc_Q20.pl fastq position\n
Author: Yang Li\n
Please feel free to contact me\n";
	exit;
	}
	
open IN, "$ARGV[0]";
open OUT, ">$ARGV[0].Q20";

my $fq;
my $n = 1;
my $i;
my $qual_last;
my $seq_len;
while(<IN>){
  chomp;
	if ($n == 1) {
		$n = 2;
		$fq .= "$_\t";}
		elsif ($n == 2) {
			$seq_len = length($_); 
			$n = 3;
			$fq .= "$_\t";}
			elsif ($n == 3) {
			$n = 4;
			$fq .= "$_\t";}
				elsif ($n == 4) {
				$fq .= "$_";
				$n = 1;
					my @seq = split /\t/, $fq;
					#print "@seq\n"; #序列变成4列信息，方便处理
					$qual_last = substr ($seq[3], -$ARGV[1]);
					my @qual = split //, $qual_last;
					#my @qual = split //, $seq[3];
					#pritn 
					$i = 0; #initial $i; 空的输出
					for (@qual) {
						my $pred_Q = ord ($_) - 33;  #将ASCII转换成Q值
						#print "pred_Q $pred_Q\n";
						if ($pred_Q < 20) {
							#print "1\n";
							#print "i is $i\n";
								my $position_cutoff = $seq_len - $i;
								substr ($seq[3], $position_cutoff) =~ s/.*//;
								substr ($seq[1], $position_cutoff) =~ s/.*//;
								last;
							}
						$i++;					    
					}
					my $seq_afterqc = join "\n", @seq;
					print OUT "$seq_afterqc\n";
 					$fq = "";
		} 
			
}

close IN;
close OUT;
