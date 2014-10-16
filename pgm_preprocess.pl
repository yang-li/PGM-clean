#!/usr/bin/perl -w

#==========================================================================
#   Author: Yang Li, Chinese Centre for Disease Control and Prevention, Beijing, China
#
#   	    <liyang@ivdc.chinacdc.cn>
#
#   Usage:
#      perl pgm_preprocess.pl [options]
#
#      Try 'perl pgm_preprocess.pl -h' for more information.
#   
#   Purpose:
#   
#   For Ion Torrent data preprocessing:
#	1, min_len, max_len, range_len
#	2, min_qual_mean, min_qual_score
#	3, ns_max_n
#   Workflow：
#	1, file format check. (fasta, fastq)
#	2, command check
#	3, execute command
#============================================================================
use strict;

use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum min max);

my $man = 0;
my $help = 0;
my %params = ('help' => \$help,
              'h' => \$help,
              'man' => \$man);
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'fastq=s',
            'fasta=s',
            'min_len=i',
            'max_len=i',
            'range_len=s',
            'min_qual_score=i',
            'min_qual_mean=i',
            'ns_max_n=i',
            'trim_left=i',
            'trim_right=i',
            'out=s',
            'phred64'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
#
##
###format check
##
#
my $file;
if(exists $params{fasta} && exists $params{fastq}) {
    &printError('fasta and fastq cannot be used together');
} elsif(exists $params{fasta} && (-e $params{fasta})) {
        #check for file format
        $file = $params{fasta};
        my $format = &checkInputFormat($file);
        unless($format eq 'fasta') {
            &printError('input file for -fasta is in '.uc($format).' format not in FASTA format');
            } 
} elsif(exists $params{fastq} && (-e $params{fastq})) {
        $file = $params{fastq};
        my $format = &checkInputFormat($file);
        unless($format eq 'fastq') {
            &printError('input data for -fastq is in '.uc($format).' format not in FASTQ format');
        }
} else {
    &printError("you did not specify an input file containing the query sequences");
}
#
##
###check command
##
#
unless( exists $params{min_len} ||
        exists $params{max_len} ||
        exists $params{range_len} ||
        exists $params{min_qual_score} ||
        exists $params{min_qual_mean} ||
        exists $params{ns_max_n} ||
        exists $params{trim_left} ||
        exists $params{trim_right} 
        ) {
    &printError('nothing to do with input data');
}
unless ( exists $params{out} ) {
    &printError('please specify a filename for output');
}

open OUT, ">$params{out}" or die "Cannot open file : $!\n";
#
##
### Core part
##
#

#order of processing:
#ns_max_n, trim_left, trim_right, min_len, max_len, range_len, min_qual_score, min_qual_mean

my $filename = $file;

while($filename =~ /[\w\d]+\.[\w\d]+$/) {
    $filename =~ s/\.([\w\d]+$)//;
    last if($filename =~ /\/[^\.]+$/);
}
open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $file |") or &printError("Could not open file $file: $!");
#data process due to file format
my ($count, $length, $seq, $seqid, $qual);
#flag
#
$count = 0;
$seq = '';
$seqid = '';
#information for data
my ($total_lines, $total_base, $num_seqs, $cur_num_seq);
my ($goodcount, $badcount, %filtercount);
if (exists $params{fasta}) {
    $total_lines =`wc -l $params{fasta}`;
    $total_lines =~ s/\s(.*)//sg;
    $num_seqs = $total_lines / 2;
    #add the funciton to display the real time progress
    while(<FILE>) {
        chomp();
        $qual='';
        if(/^>(.*)/o) {
            $length = length($seq);
            if($length) {
                &processEntry($length,$seq,$seqid,$qual);
            }
            $seqid = $1;
            $seq = '';
        } else {
                $seq .= $_;
            }
        }
    &processEntry($length,$seq,$seqid,$qual);
}elsif (exists $params{fastq}) {
    $total_lines =`wc -l $params{fasta}`;
    $total_lines =~ s/\s(.*)//sg;
    $num_seqs = $total_lines / 2;
    while(<FILE>) {
        chomp();
        if($count == 0 && /^\@(.*)/o) {
            $length = length($seq);
            if($length) {
                &processEntry($length,$seq,$seqid,$qual);
            }
            $seqid = $1;
            $seq = '';
            $qual = '';
            } elsif($count == 1) {
            $seq = $_;
            } elsif($count == 3) {
            $qual = $_;
            $count = -1;
            }
            $count++;
        }
       &processEntry($length,$seq,$seqid,$qual);
    }

print "The information for $file:\n";
printf "%20s%10d\n",'Number of reads:',"$num_seqs";
printf "%20s%10d\n",'Number of bases:',"$total_base";
print "Process Done\nReads filted by the following commands:\n";
printf "%10s%10s%10s\n",'Command','Params','NumReads';
while (my ($key, $value) = each %filtercount) {
    printf "%10s%10s%10s\n","$key","$params{$key}","$value";
}

#==========================subs=================================
sub checkInputFormat {
    my $temp_file = shift;
    my ($format,$count,$fasta,$fastq);
    $count = 3;
    $fasta = $fastq  = 0;
    $format = 'unknown';

#    open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $temp_file |") or &printError("Could not open file $file: $!");
    open(FILE,"perl -pe 's/\r\n|\n/\n/g' < $temp_file |") or &printError("Could not open file $file: $!");
    while (<FILE>) {
#        chomp();
 #       next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/o) {
            $fasta = 1;
        } elsif($fasta == 1 && (/^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o)) {
            $fasta = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/o) {
            $fastq = 1;
        } elsif($fastq == 1 && (/^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o)) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+/o) {
            $fastq = 3;
        }
    }
    close(FILE);
    if($fasta == 2) {
        $format = 'fasta';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }
    return $format;
}

sub checkRange {
    my ($range,$val) = @_;
    my @ranges = split(/\,/,$range);
    foreach my $r (@ranges) {
	my @tmp = split(/\-/,$r);
	return 0 if($val < $tmp[0] || $val > $tmp[1]);
    }
    return 1;
}

sub processEntry {
    my ($length,$seq,$seqid,$qual) = @_;
    #check that sequence and quality are same length
    if(defined $qual && length($qual) && $length != length($qual)) {
        &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"");
    }
    #remove anything non-alphabetic from sequences
    $seq =~ tr/a-zA-Z//cd;
    $total_base += $length;
    &processData($seqid,$seq,$qual);
}


sub printError {
    my $msg = shift;
    print STDERR "\nERROR: ".$msg.".\n\nTry \'perl pgm_preprocess.pl -h\' for more information.\nExit program.\n";
    exit(0);
}

sub processData {
    my ($sid,$seq,$qual) = @_;
    #assume sequence is good ;-)
    my $good = 1;
    my $seqn = uc($seq);
    my $qualn = $qual;
    my $begin = 0;
    my $end = 0;
    my ($length,$bylength,$qualsnums);
    
    #check for N's in sequence
    if($good && exists $params{ns_max_n}) {
        my $ns = ($seqn =~ tr/N//);
        if($good && exists $params{ns_max_n} && $ns > $params{ns_max_n}) {
            $good = 0;
            $filtercount{ns_max_n}++;
        }
    }

    #trim sequence ends
    if($good && exists $params{trim_left}) {
        $begin += $params{trim_left};
        if($begin >= length($seqn)) {
            $good = 0;
            $filtercount{trim_left}++;
        } else {
            $seqn = substr($seqn,$begin);
            $qualn = substr($qualn,$begin) if(defined $qualn && length($qualn));
        }
    }
    if($good && exists $params{trim_right}) {
        $end += $params{trim_right};
        $length = length($seqn);
        if($end >= $length) {
            $good = 0;
            $filtercount{trim_right}++;
        } else {
            $seqn = substr($seqn,0,$length-$end);
            $qualn = substr($qualn,0,$length-$end) if(defined $qualn && length($qualn));
        }
    }

    #check if trim to certain length
    $length = length($seqn);
    if($good && exists $params{trim_to_len} && $length > $params{trim_to_len}) {
        $seqn = substr($seqn,0,$params{trim_to_len});
        $qualn = substr($qualn,0,$params{trim_to_len}) if(defined $qualn && length($qualn));
        $end += ($length-$params{trim_to_len});
    }

    #check for sequence length
    $length = length($seqn);
    $bylength = ($length ? 100/$length : 0);
    if($bylength == 0) {
        $good = 0;
        $filtercount{zero_length}++;
    }
    if($good && exists $params{min_len} && $length < $params{min_len}) {
        $good = 0;
        $filtercount{min_len}++;
    }
    if($good && exists $params{max_len} && $length > $params{max_len}) {
        $good = 0;
        $filtercount{max_len}++;
    }
    if($good && exists $params{range_len} && !&checkRange($params{range_len},$length)) {
        $good = 0;
        $filtercount{range_len}++;
    }

    #check for quality scores
    if($good && defined $qualn && (exists $params{min_qual_score} || exists $params{min_qual_mean} )) {
        my ($err);
        if($qualsnums) {
            if($begin > 0) {
                shift(@$qualsnums) foreach(1..$begin);
            }
            if($end > 0) {
                pop(@$qualsnums) foreach(1..$end);
            }
        } else {
            if(exists $params{phred64}) { #scale data to Phred scale if necessary
                ($qualsnums,$err) = &convertQualAsciiToNumsPhred64($qualn);
                if($err) {
                    &printError("The sequence quality scores are not in Phred+64 format");
                }
            } else {
                $qualsnums = &convertQualAsciiToNums($qualn);
            }
        }
        if($good && exists $params{min_qual_score} && min(@$qualsnums) < $params{min_qual_score}) {
            $good = 0;
            $filtercount{min_qual_score}++;
        }
        if($good && exists $params{min_qual_mean} && &getArrayMean(@$qualsnums) < $params{min_qual_mean}) {
            $good = 0;
            $filtercount{min_qual_mean}++;
        }
    }
    #print OUT
        if ($good && (exists $params{fasta})) {
            print OUT ">$seqid\n$seq\n";    
        } elsif ($good && (exists $params{fastq})) {
            print OUT "@$seqid\n$seq\n+\n$qualn\n";
            }   
}

sub convertQualAsciiToNums {
    my $qual = shift;
    my @nums;
    my @ascii = split(//,$qual);
    foreach(@ascii) {
	push(@nums,(ord($_) - 33));
    }
    return \@nums;
}

sub convertQualAsciiToNumsPhred64 {
    my $qual = shift;
    my @nums;
    my $tmp;
    my $err = 0;
    my @ascii = split('',$qual);
    foreach(@ascii) {
        $tmp = (ord($_) - 64);
        if($tmp < 0) {
	    $err = 1;
	    last;
	}
	push(@nums,$tmp);
    }
    return (\@nums,$err);
}


sub getArrayMean {
    return @_ ? sum(@_) / @_ : 0;
}

__END__
 
=head1 NAME
 
pgm_preprocess - data preprocessing for ion torrent Platform
 
=head1 SYNOPSIS
 
perl pgm_preprocess.pl [-h] [-help] [-man] [-fastq input_fastq_file] [-fasta input_fasta_file] [-min_len int_value] [-max_len int_value] [-range_len ranges] [-min_qual_score int_value] [-min_qual_mean int_value] [file ...] [-trim_left int_value] [-trim_right int_value] [-phred64] [-out filename]
 
 Options:
   -help|h            brief help message
   -man               full documentation
 
=head1 OPTIONS
 
=over 8
 
=item B<-help> | B<-h>
 
Print a brief help message and exits.
 
=item B<-man>
 
Prints the manual page and exits.

=item B<***** INPUT OPTIONS *****>

=item B<-fastq> <file>

Input file in FASTQ format that contains the sequence and quality data.

=item B<-fasta> <file>

Input file in FASTA format that contains the sequence data. 

=item B<-phred64>

Quality data in FASTQ file is in Phred+64 format (http://en.wikipedia.org/wiki/FASTQ_format#Encoding). Not required for Illumina 1.8+, Sanger, Roche/454, Ion Torrent, PacBio data.

=item B<***** OUTPUT OPTIONS *****>

=item B<-out> <string>

By default, the output files are created in the same directory as the input file containing the sequence data with the name defined by user. To change the output filename and location, specify the filename using this option. 

=item B<***** FILTER OPTIONS *****>

=item B<-min_len> <integer>

Filter sequence shorter than min_len.

=item B<-max_len> <integer>

Filter sequence longer than max_len.

=item B<-range_len> <string>

Filter sequence by length range. Multiple range values should be separated by comma without spaces.

Example: -range_len 30-70,150-250

=item B<-min_qual_score> <integer>

Filter sequence with at least one quality score below min_qual_score. (recommend: 10)

=item B<-min_qual_mean> <integer>

Filter sequence with quality score mean below min_qual_mean. (recommend: 17)

=item B<-ns_max_n> <integer>

Filter sequence with more than ns_max_n Ns. (recommend: 1)

=item B<***** TRIM OPTIONS *****>

=item B<-trim_left> <integer>

Trim sequence at the 5'-end by trim_left positions.

=item B<-trim_right> <integer>

Trim sequence at the 3'-end by trim_right positions.

=item B<***** ORDER OF PROCESSING *****>

ns_max_n, trim_left, trim_right, min_len, max_len, range_len, min_qual_score, min_qual_mean

=item B<->
 
=back
 
=head1 AUTHOR

Yang Li, C<< <liyang@ivdc.chinacdc.cn> >>
 
=cut

