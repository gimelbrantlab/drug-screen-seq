#! perl  
####################################################################
# v.7.2	12/20/18	Sasha - fixing issues with UMI calling for Clara (HONDA sequence run)
# v.7.1	7/20/18	Sasha - corrected an issue with incorrect detection of asymmetric flag
# v.7		6/30/17	Sasha - integration of UMI and SNP analysis
#					   yaml config updated to v.6.0
#					   TryRevComp removed
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
#### input: SAM file for specific gene region/target
#### IMPORTANT: sorted by the read# (so that the R1 and R2 reads are next to each other)
#### output: table of SNPs per 
####################################################################

use YAML::XS 'LoadFile';
use feature ":5.10";
no warnings 'experimental';

my $REVERSE = 0x10;        ## in FLAG, 5th bit is set if revcomp
my $SUPPLEMENTARY = 0x800;	## in FLAG, this is set for added short read in as2.pl

### load params from command line ###############################
$yaml=shift; ### yaml file
unless ($yaml ne ''){die "$0: <yaml_file> is needed\n"};  ## TODO: better error checking for YAML loading
my $config = LoadFile($yaml);
## TODO: check success of yaml loading

my $infile =shift;  ## sorted SAM file for this gene
my $gene_name=shift;
my $timestamp=shift;

# access general yaml content ###########################################
$ExpectedConfigVersion = "6.0";
$ConfigVersion=$config->{version};
unless($ExpectedConfigVersion eq $ConfigVersion) {die "$0: YAML config file version $ConfigVersion. Expecting $ExpectedConfigVersion\n"};

$DEBUG = $config->{diagnostics}->{DEBUG};
if ($config->{diagnostics}->{OnlyOutputDiagnostics} eq "yes") {$OnlyOutputDiagnostics=1}
else {$OnlyOutputDiagnostics=0};

$MAPQmin=$config->{step2}->{minimal_MAPQ}; # min quality -- global
$LongReadSplitPattern=$config->{step2}->{LongReadSplitPattern};
$ShortReadSplitPattern=$config->{step2}->{ShortReadSplitPattern};

### TODO: the following will not be needed ~~~~~~~~~~~~~~~~~
$globalTLENmin=$config->{step2}->{insert_len_min};
$globalTLENmax=$config->{step2}->{insert_len_max};
$amplicon_size_tolerance=$config->{step2}->{amplicon_size_tolerance};
$SNPnoUMIpattern=$config->{step2}->{SNPnoUMIpattern};
$T2_small=$config->{step2}->{T2_small};
$T2_small_rc = &revcomp($T2_small);
#~ if ($config->{diagnostics}->{TryRevComp} eq "yes") {$TryRevcomp=1}
#~ else {$TryRevcomp=0};

$flankLeft	= $config->{genes}->{$gene_name}->{SNPFlank_left};
#### This is temporary: !!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		$snp = $flankLeft;
		#~ $snp_revcomp = revcomp($snp);
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## access gene-specific YAML content #####################################
$SNP	=$config->{genes}->{$gene_name}->{SNP};
	if($SNP =~ m/(\w+)\[\w+\/\w+\](\w+)/gi)
		{$SNPleftFlank=$1; $SNPrightFlank=$2;}
	else {die "\nError in $0:\n problem with SNP parsing >$SNP<\n"};
	
	#~ die "\n$SNP\n$SNPleftFlank $SNPrightFlank\n";

## The following doesn't matter for the short read (R1)
#~ $size=$config->{genes}->{$gene_name}->{amplicon_size};
#~ $TLENmin=$size - $amplicon_size_tolerance; 
   #~ if($TLENmin < $globalTLENmin){$TLENmin=$globalTLENmin};
#~ $TLENmax=$size + $amplicon_size_tolerance;
   #~ if($TLENmax > $globalTLENmax){$TLENmax=$globalTLENmax};
   

# load samples_barcodes ###############################################
for (sort keys %{$config->{samples_barcodes}}) {
	$descr=$_; $bc=$config->{samples_barcodes}->{$_};
	push (@samples_barcodes,[$descr, $bc]);
	$expected_barcodes{$bc}++; ## for counting unexpected bcodes
}

#### Now deal with the SAM file  #########################################
open(S,"<$infile") or die("$0: cannot open $infile\n");
### init counters
$read_pairs=
$MAPQ_low=
$T2_not_found=
$bc_mismatch=
$accum=
$extraUMIaccum=
$QNAME_mismatch=
$unknown_bc=
$flank_not_found=0;

while(<S>){
	if(/^@/){next};
	### Here, line1 and line2 are just the order of reading; no assumptions are made which one had been added in asymm process
	($line1{QNAME}, $line1{FLAG}, $line1{RNAME}, $line1{POS}, $line1{MAPQ}, $line1{CIGAR}, $line1{RNEXT}, $line1{PNEXT}, $line1{TLEN}, $line1{SEQ}, $line1{QUAL}) = split;
	$_=<S>;	
	($line2{QNAME}, $line2{FLAG}, $line2{RNAME}, $line2{POS}, $line2{MAPQ}, $line2{CIGAR}, $line2{RNEXT}, $line2{PNEXT}, $line2{TLEN}, $line2{SEQ}, $line2{QUAL}) = split;

	
## IDs should be identical for r1 and r2
	if($r1{QNAME} ne $r2{QNAME}) {$QNAME_mismatch++; $_=<S>; next}; 
	##  very crappy.  TODO: change R1/R2 matching correction so that it doesn't skip data
	$read_pairs++;

### ~~~ which one is short read
	if	(($line1{FLAG}+0) & $SUPPLEMENTARY){  ## line1 is the short read added in asymmetric process
		%R1=%line1; %R2=%line2;
		}
	elsif (($line2{FLAG}+0) & $SUPPLEMENTARY){  ## line2 is the short read added in asymmetric process
		%R1=%line2; %R2=%line1;
		}
	else{ 
		print STDERR "$0: No asymm flag found: expected FLAG=$SUPPLEMENTARY   R1=$r1{FLAG}   R2=$r2{FLAG}\n";
		next
		};

	if(($R2{MAPQ} < $MAPQmin) ){$MAPQ_low++; next};  
	
## reads mapped in revcomp orientation?  set the flag
	#~ if (($R2{FLAG}+0) & $REVERSE) {$R2{SEQ} = revcomp($R2{SEQ});}
	if (($R2{FLAG}+0) & $REVERSE) {$revR2=1} else{$revR2=0};
	$revR1=0;
	#::: R1 is never reverse in asymmetric reads::: if (($R1{FLAG}+0) & $REVERSE) {$revR1=1} else{$revR1=0}   
	if($revR2){$R2{SEQ} = &revcomp($R2{SEQ})}
	
### now see what the barcodes, UMIs, and SNPs are #########################
    ## short read (with UMI)
	($R1{spacer}, $R1{bc}, $R1{tail}, $R1{UMI}, $R1{rest}) = unpack($ShortReadSplitPattern, $R1{SEQ});
  
   ## long read
	($R2{spacer}, $R2{bc}, $R2{tail}, $R2{rest}) = unpack($LongReadSplitPattern,  $R2{SEQ});
   
## Identify what the barcode set is
	#~ my $theBarCodePair=&getTheBarcode($R1{bc}, $revR1, $R1{tail}, $R2{bc},$revR2, $R2{tail}, $T2_small, $T2_small_rc);
	my $theBarCodePair=&getTheBarcode($R1{bc}, $R1{tail}, $R2{bc}, $R2{tail}, $T2_small);


	if($theBarCodePair eq "none"){$T2_not_found++; next};
	unless(exists($expected_barcodes{$theBarCodePair})){$unknown_bc++; next};

	$all_processed++; ### This is counting all assessed read pairs, NOT including those with "undeclared" barcodes

## Skip if UMI already encountered (for given barcode and gene)	
	if(exists($data{$theBarCodePair}{$R1{UMI}}))
		{$data{$theBarCodePair}{extraUMIs}++;
		next;
		}
		else {$data{$theBarCodePair}{$R1{UMI}}=1;  
		};

### now do the "genotyping" (only if unique UMI) ###############################

	if(	$R2{rest}=~m/$SNPleftFlank(.)$SNPrightFlank/gi)
		{$gt=uc($1)}
	else	{$gt="X"};	
	
									#~ $FOO++; if($FOO > 50){ exit; }
									#~ print "$theBarCodePair	Flag: $R2{FLAG}	R2: $R2{rest}	SNP: $SNPleftFlank.$SNPrightFlank	call=$gt	$R1{QNAME}\n";
	given ($gt) {
		when ("A") {$data{$theBarCodePair}{A}++}
		when ("C") {$data{$theBarCodePair}{C}++}
		when ("G") {$data{$theBarCodePair}{G}++}
		when ("T") {$data{$theBarCodePair}{T}++}
		when ("N") {$data{$theBarCodePair}{N}++}
		when ("X") {$data{$theBarCodePair}{X}++;  ### NO SNP flank found
					$flank_not_found++;
					#~ print STDERR "$0";
					## TODO reset UMI found flag:  maybe the next read with this UMI will have correct SNP flanks
					}  
		default {next};
	};
}
# we are done reading data and counting, now let's print some output


################################################
# output results -- all on one line per barcode; the following fields:
#~ 1  	timestamp
#~ 2 	gene_id
#~ 3	sample_id
#~ 4	count_a
#~ 5	count_t
#~ 6	count_g
#~ 7	count_c
#~ 8	count_n
#~  9	count_no_flank
#~ 10	log
################################################
{my $i;
	$bc_mismatch=0;
# loop through and print output
	for($i=0; $i<=($#samples_barcodes); $i++){
		$output_line="$timestamp\t$gene_name\t$samples_barcodes[$i][0]";

	if (exists($data{$samples_barcodes[$i][1]})){
	       ## count
		unless ($OnlyOutputDiagnostics){
			if (exists($data{$samples_barcodes[$i][1]}{A})){$tmp=sprintf ("\t$data{$samples_barcodes[$i][1]}{A}");$output_line.=$tmp} else {$output_line.="\t0"};
			if (exists($data{$samples_barcodes[$i][1]}{T})){$tmp=sprintf ("\t$data{$samples_barcodes[$i][1]}{T}");$output_line.=$tmp} else {$output_line.="\t0"};
			if (exists($data{$samples_barcodes[$i][1]}{G})){$tmp=sprintf ("\t$data{$samples_barcodes[$i][1]}{G}");$output_line.=$tmp} else {$output_line.="\t0"};
			if (exists($data{$samples_barcodes[$i][1]}{C})){$tmp=sprintf ("\t$data{$samples_barcodes[$i][1]}{C}");$output_line.=$tmp} else {$output_line.="\t0"};
			if (exists($data{$samples_barcodes[$i][1]}{N})){$tmp=sprintf ("\t$data{$samples_barcodes[$i][1]}{N}");$output_line.=$tmp} else {$output_line.="\t0"};
			if (exists($data{$samples_barcodes[$i][1]}{X})){$tmp=sprintf ("\t\{$data{$samples_barcodes[$i][1]}{X}\}");$output_line.=$tmp} else {$output_line.="\t{0}"};
			if (exists($data{$samples_barcodes[$i][1]}{extraUMIs})){$tmp=sprintf ("\t<$data{$samples_barcodes[$i][1]}{extraUMIs}>");$output_line.=$tmp} else {$output_line.="\t<0>"};
			
		};       
		$accum +=	$data{$samples_barcodes[$i][1]}{A}+
					$data{$samples_barcodes[$i][1]}{T}+
					$data{$samples_barcodes[$i][1]}{C}+
					$data{$samples_barcodes[$i][1]}{G};
					
		$extraUMIaccum += 	$data{$samples_barcodes[$i][1]}{extraUMIs};		
	       }
	else { ### no data found
		unless ($OnlyOutputDiagnostics){$output_line.="\t0\t0\t0\t0\t0\t{NA}\t<NA>";}
		$bc_mismatch++;}

	### now print the assembled line
	print "$output_line\n";
	}
}	
##################

if($DEBUG){
	################################################
	if($all_genotyped == 0){$all_genotyped=-1};
	if($read_pairs == 0){$read_pairs=-1};
	
	print "#\n#	~~~~~~~~~ Debug info. To disable this output, set diagnostics:DEBUG to 0\n";
	print "#	From	$infile\n#	$read_pairs	read pairs\n";
			$sb=$#samples_barcodes+1;	
	print "#	$sb	barcodes declared; of these, $bc_mismatch had no reads found\n";

	print "#	~~~~~~~~~ Filtered out:\n";
	print "#	$QNAME_mismatch	unpaired reads\n";
	print "#	$MAPQ_low	MAPQ<$MAPQmin\n";
	#~ print  "#	$TLEN_same_sign\treads aligned same direction\n";
	print "#	$T2_not_found	T2_not_found\n";
				$tmp=sprintf("%.1f", 100*$unknown_bc/$read_pairs); 
	print "#	$unknown_bc	reads with undeclared barcodes ($tmp% of $read_pairs)\n";			
	
	print "#	~~~~~~~~~ Processed:\n";
	if($all_processed == 0) {$all_processed=-1};
	if($accum == 0) {$accum=-1};
	
	
			$tmp=sprintf("%.1f", 100*$all_processed/$read_pairs);	
	print "#	$all_processed	Read pairs that passed all the above filters ($tmp% of $read_pairs)\n";
	
			$sb1=$sb-$bc_mismatch; 		
			$tmp=sprintf("%.1f", 100*$accum/$all_processed); 
	print "#	$accum	<<<<< Total unique genotype calls ($tmp% of $all_processed) among $sb1 of declared barcodes with reads\n";
	
			$tmp=sprintf("%.1f", 100*$flank_not_found/$accum);
	print "#	{$flank_not_found}	SNP flank not found ($tmp% of $accum)	~~~ $SNPleftFlank [_] $SNPrightFlank\n";

			$tmp=sprintf("%.1f", 100*$extraUMIaccum/$all_processed);
	print "#	<$extraUMIaccum>	extra UMIs ($tmp% of $all_processed)\n";
	
	print "\n";

	#~ if($TryRevcomp){
		#~ if($all_genotypedRC == 0){$all_genotypedRC=-1};
		#~ if($read_pairsRC == 0){$read_pairsRC=-1};

		#~ print "~~~~~~~ RevComp for SNP flank\n";	
		#~ $tmp=sprintf("%.1f", 100*$all_genotypedRC/$read_pairs);	
		#~ print "$all_genotypedRC	Read pairs that passed all the above filters ($tmp% of total)\n";
		
		#~ $tmp=sprintf("%.1f", 100*$accumRC/$all_genotypedRC); 
		#~ print  "$accumRC	Total unique genotype calls ($tmp% of passed) for declared barcodes with reads\n";
	#~ }
}

###############################################################
###############################################################
###############################################################

sub revcomp {      #from http://www.perlmonks.org/?node_id=197793
  my $dna = shift; 
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}


####### updated 12/21/18
####### we know that the seq was assymmetrical; R1 is short read and R2 is long read
sub getTheBarcode{ 
	my($R1_bc, $R1_tail, $R2_bc, $R2_tail, $T2)=@_;
	my $answer="NA";
## sanity check: R1 (short read) should have T2 ####################################
	if($R1_tail =~ m/$T2/) { 	$answer = "$R2_bc"."_"."$R1_bc" ;} else {$answer="none";}

	return $answer;	
}


	
__END__
