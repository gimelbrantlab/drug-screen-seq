#! perl  
####################################################################
# v.0.2	2/22/17	Sasha - deal with asymmetric reads
####################################################################
# v.6.2  	1/4/17	Sasha - more diagnostics in DEBUG
# v.6  	12/30/16	Sasha - deal with the short R1 read (barcode only; not alignable to genome)
# v.5.1	12/26/16	Sasha - unified YAML5.0 with UMI branch
# v.5.0	12/3/16	Sasha - new output format, per line
# v.0.10	9/30/16	Sasha - YAML4.5; additional spacer
# v.0.9	8/26/16	Sasha - multiple updates to YAML (v.4.0)
# v.0.8	8/16/16	Sasha - use revcomp (YAML v.3.1)
# v.0.7	8/2/16	Sasha - SNP counting
# v.0.6	6/24/16	Sasha - clean up code; output libraries in sorted order
# v.0.5	3/24/16	Sasha - use YAML
# v.0.4	3/10/16	Sasha - cleaned up code a bit
# v.0.3	3/7/2016	Sasha - arbitrary barcode pairs
# v.0.2	3/6/2016	Sasha - make more parameters
# v.0.1 	3/4/2016	Sasha - initial version
####################################################################
#### input: SAM file (from long fastq) in the same order as fastq; FASTQ for short reads 
#### output: only read pairs where the long read aligned
####################################################################

#~ use YAML::XS 'LoadFile';
use feature ":5.10";
no warnings 'experimental';

### const
$EOF=-65;
my $mask = 1 << 4;  ## in FLAG, 5th bit is set if revcomp

### load params from command line ###############################
my $SAMfile =shift;  ## SAM file 
my $FQfile = shift; ## short FASTQ

#### Now deal with the input files  #########################################
open(S,"<$SAMfile") or die("$0: cannot open $SAMfile\n");
open(F,"<$FQfile") or die("$0: cannot open $FQfile\n");

### init counters
$read_pairs=$QNAME_mismatch=$accum=0;

while(<S>){
	if (/^@/) {print; next}; ## skip @ lines in SAM file
#### r2 -- from the sam file
	($r2{QNAME}, $r2{FLAG}, $r2{RNAME}, $r2{POS}, $r2{MAPQ}, $r2{CIGAR}, $r2{RNEXT}, $r2{PNEXT}, $r2{TLEN}, $r2{SEQ}, $r2{QUAL}) = split;
	$line2=$_;
	#~ die (">>>from sam=$r2{QNAME}\n$line2\n");

#### r1 -- from the fastq file	
	$ans=readFQ(); ## read from fastq file	
	if($ans eq $EOF){exit(0)};
	if($ans eq -1){die "problem with FASTQ reading; returned $ans\n"};
	
	$read_pairs++;

## IDs should be identical for R1 and R2
	if($r1{QNAME} ne $r2{QNAME}) {$QNAME_mismatch++; $_=<S>; print STDERR "mismatch: $read_pairs lines \n$r1{QNAME}\n$r2{QNAME}\n";next}; ## <-- very crappy.  TODO: change R1/R2 matching correction so that it doesn't skip data

## skip pairs whete none of the reads aligned
	if(($r2{RNAME} eq "*")){$no_aln++; next };  # not aligned
### ~~~ DEAL with short R1
	$foo=length($r1{SEQ});
	$r1{CIGAR}=sprintf("%iM",$foo);
	
	### set flag by toggling revcomp flag in r2
	#~ my $mask = 1 << 4;  ## in FLAG, 5th bit is set if revcomp
	
	$r1{FLAG}=(($r2{FLAG}+0) ^ $mask)+2048;  ## +0 to force conversion to number 
	
	$line1=sprintf("$r1{QNAME}	$r1{FLAG}	$r2{RNAME}	$r2{POS}	$r2{MAPQ}	$r1{CIGAR}	$r2{RNEXT}	$r2{PNEXT}	$r2{TLEN}	$r1{SEQ}	*");

	print "$line2"; print "$line1\n";   ## long then short
}

sub readFQ{
	my $foo=<F>;
	unless($foo){return $EOF};
	unless($foo=~/^@/){ return -1};
	$foo=~s/@//; ## remove the "@" from the name
	($r1{QNAME})=split(' ',$foo);
		#~ $r1{QNAME}=$foo;
	$r1{SEQ}=<F>; chomp $r1{SEQ};
	$foo=<F>; ## skip the next two lines
	$foo=<F>;
	unless($foo){return $EOF};
	#~ die "$r1{QNAME} ~~~ $r1{SEQ}\n";
	#~ $r1{QNAME}="FOUND";
	return 0;
}

__END__

