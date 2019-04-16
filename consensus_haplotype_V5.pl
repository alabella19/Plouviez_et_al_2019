#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
#use Bio::SeqIO;
#use Bio::SimpleAlign;
#use Bio::Tools::Run::StandAloneBlast;
#use Bio::Tools::Run::StandAloneBlastPlus::BlastMethods;
#use Bio::Tools::dpAlign; 
use List::Util 'max';
use List::Util qw(first);

##This script is designed to take a process_____.nex file and combine it with the original ______.fas file and the Individual

###Now we must account for the gaps between se

##VERSION 4: 9.30.14 Fixes some weird problem with the creation of sequences. It started to guess that consensus sequences were PROT not DNA


usage() if(@ARGV <2 or ! GetOptions('proc_file=s'=> \$proc_file, 'fas_file=s'=> \$fas_file, 'indv=s'=>\$indv, 'nex_file=s'=>\$nex));
sub usage{
	print "\n USAGE: Must have process file from R-script and the fasta file that was input to the SNP caller\n -proc_file <process____.nex>\n -fas_file <_____.fas>\n";
	exit;
}

##Parse the NEX file to get all of the SNPs
open (INPUT, $nex) || die $!;
@nex_array = <INPUT>;
close(INPUT);

$SNPS = $nex_array[6];

$SNPS =~ s/\t//g;
$SNPS =~ s/\n//;
$SNPS =~ s/\,//g;
@SNPS = split('\s',$SNPS);
@SNPS = sort {$a <=> $b} @SNPS;
#print "SNPS Starts jere @SNPS\n";

#print Dumper (\@SNPS);
#exit;

###Parse the process file to get the SNP positions and the haplotypes and the reads!

open (INPUT, $proc_file) || die $!;
@proc_array = <INPUT>;
close(INPUT);

$has_three = 0;
if($proc_file =~ /^3/){
	$has_three = 1;
}

$num_block = $proc_file;
$num_block=chop($num_block);

$indv=~s/\///;
$indv=~s/\-/\./;

$outfile=$fas_file;
$outfile=~s/\.fasta//;
$outfile=~s/\.fas//;
$outfile=$outfile.".haps."."$num_block.fasta";
if($has_three == 1){
	$outfile="3".$outfile;
}

$line1 = shift(@proc_array);

$line1 =~ s/\n//;
$line1 =~ s/\r//;
$line1 =~ s/\.//g;

@line_array = split('\s+',$line1);
#print Dumper(\@line_array);
$line_size = scalar(@line_array);
#print "len : $line_size";
@pos = ();
$n = 2;
while($n < $line_size){
	
		push(@pos, $line_array[$n]);
		$n = $n +2;
}

#print Dumper (\@pos);

%hap_hash1 = ();
%hap_hash2 = ();
%hap_hash3=();

$line2= shift(@proc_array);
$line2 =~ s/\n//;
$line2 =~ s/\r//;

@snps = split('\s+',$line2);
$hap1_name = shift(@snps);
##FOR test3rd
$hapT1 = join('',@snps);
$hap1_name =~ s/haplotype/$indv\./;
#print Dumper (\@snps);
#print "SNP:-$snps[1]-";
$n = 0; 
foreach $char (@snps){
	$hap_hash1{$pos[$n]} = $char;
	$n = $n +1;
}

$line3= shift(@proc_array);
$line3 =~ s/\n//;
$line3 =~ s/\r//;
@snps = split('\s+',$line3);
$hap2_name = shift(@snps);
##FOR test3rd
$hapT2 = join('',@snps);
$hap2_name =~ s/haplotype/$indv\./;
$n = 0; 
foreach $char (@snps){
	$hap_hash2{$pos[$n]} = $char;
	$n = $n +1;
}
if($has_three == 1){
	$line4 = shift(@proc_array);
	$line4 =~ s/\n//;
	$line4 =~ s/\r//;
	@snps = split('\s+',$line4);
	$hap3_name = shift(@snps);
	##FOR test3rd
	$hapT3 = join('',@snps);
	$hap3_name =~ s/haplotype/$indv\./;
	$n = 0; 
	foreach $char (@snps){
		$hap_hash3{$pos[$n]} = $char;
		$n = $n +1;
	}
}

#print "$hapT1\n$hapT2\n$hapT3\n";
if($has_three == 1){
	test3rd();
}
if($RECOMB != 0){
	$has_three = 0;
	#get rid of recombinant allele
	if($RECOMB == 1) {
		%hap_hash1 = %hap_hash3;
		$hap1_name = $hap3_name;
	}
	if($RECOMB == 2) {
		%hap_hash2 = %hap_hash3;
		$hap2_name = $hap3_name;
	}
}
	
	
#print "Recombinant is $RECOMB";
#exit;

###Get all of the reads
%read_hash = ();
foreach $line (@proc_array){
	@read = split('\s+', $line);
	$read = @read[0];
	$read_hash{$read} = 1;
}

#print Dumper (\%read_hash);


##MAKE A SINGLE COnSENSUS FILE FROM THE FASTA

open (INPUT, $fas_file) || die $!;
@fas_array = <INPUT>;
close(INPUT);

@align_data = ();
foreach $line (@fas_array){
	if($line !~ />/){
		$n_line = $line;
		$n_line =~ s/\-/N/g;
		$n_line =~ s/\?/\-/g;
		#$n_line =~s/N/\-/g;
		push(@align_data, $n_line);
	}else{
		push(@align_data, $line);
	}
	
}

open FILE, ">", "align_data.fasta";
print FILE @align_data;
close FILE;

###get the full consensus sequence
$all_cons_fasta = Bio::AlignIO->new(-file => "align_data.fasta");
$all_cons_algn = $all_cons_fasta->next_aln();
#$all_cons_data = Bio::SimpleAlign->new($all_cons_algn);
###We were missing leading bases that were in a low frequency 
$all_cons = $all_cons_algn->consensus_iupac(10);
##Change this to get at the coverage of the non-linked regions 
$CONS = $all_cons_algn->consensus_string(10);
#print "Cons is $CONS\n";
$CONS_block = $all_cons_algn->consensus_string(1);
#CHANGED THIS TO ACCOUNT FOR Ns
$all_cons =~ s/\?/N/g; 	
#print "allCons: $all_cons\n\n";

#get just the reads that are in the current block...
$current_read = "";
@block_array = ();
foreach $line (@align_data){
	if($line =~ />/){
		$current_read = $line;
		$current_read =~ s/>//;
		$current_read =~ s/\n//;
		if($read_hash{$current_read} == 1){
			push(@block_array, ">".$current_read."\n");
			#print "yes: $current_read\n";
		}else{
			#print "no: $current_read\n";
		}
	}elsif($read_hash{$current_read} == 1){
		push(@block_array, $line); 
	}
}
open FILE, ">", "align_data.fasta";
print FILE @block_array;
close FILE;
#exit;		
$block_cons_fasta = Bio::AlignIO->new(-file => "align_data.fasta");
$block_cons_algn = $block_cons_fasta->next_aln();
#$all_cons_data = Bio::SimpleAlign->new($all_cons_algn);
###not using IUPAC for this... but maybe we should? We might lose low frequency reads... but 
##to account for low freq reads you can change the parameter from 1 to 10?
$block_cons = $block_cons_algn->consensus_string(1);
#CHANGED THIS TO ACCOUNT FOR Ns
$block_cons =~ s/\?/N/g; 
##JUST FOR TESTING
#$block_cons = substr($block_cons, 55, 100);

#print "block consensus $block_cons\n\n";

open FILE, ">", "align_data_2.fasta";
print FILE ">allConsensus\n$all_cons\n";
print FILE ">blockConsensus\n$block_cons\n";
close FILE;
#exit;
`muscle -in align_data_2.fasta -out align_data_3.fasta`;
$in = Bio::AlignIO->new(-file => 'align_data_3.fasta' , -format => 'fasta');
$align = $in->next_aln();
#$align = Bio::SimpleAlign->new(-file =>'align_data_3.fasta');
$offset_cons = $align->consensus_string(51);
#print "Consensus is: $offset_cons\n";
#exit;
###Now I have to figure out how many ???s are in front of the sequences. 
@offset = split('',$offset_cons);
#print "\n".$offset[1];
$n=0;
$end = 0;
while($end == 0){
	$test = shift(@offset);
	#print $test;
	if($test !~ /\?/ && $test !~ /N/){
		$end =1;
		last;
	}else{
		$n = $n +1;
	}
}

#print "Offset is $n\n";
$offset = $n; 


###NOW GET THE SNP POS! 

#WHAT IF THIS IS THE LAST SNP!?
#print keys %hap_hash1;
$last_SNP = max (keys %hap_hash1);
#print "Max Pos is $last_SNP\n";
$ind = scalar keys %hap_hash1;
#print "INDEX IS $ind\n";
$n_snps = scalar(@SNPS) - 2;
#print Dumper(@SNPS);
#print "L = $n_snps\n";
if($ind == $n_snps){
		#print "Is LAST\n";
		$end_SNP = length($CONS);
		$end_SNP = $end_SNP -1; ##to account for 0 start
}else{
	#$ind = $ind; # to get to the next snp
	$end_SNP = $SNPS[$ind];
	#print "this the index $ind for this hash\n";
	#print $SNPS[0]."\t".$SNPS[1];
	#print Dumper sort @SNPS;
	##SOMETHING IS WRONG WIHT SNPS IT IS MISSING 59
	$end_SNP = $end_SNP - 1; ###account for one base before the SNP
}
#print "END SNP is ".$end_SNP."\n";

####NOW get the part of the consensus sequence we want to use from $offset to $end_SNP

$len = $end_SNP - $offset; 

#print "len is $len\n";

$CONS = substr($CONS, $offset, $len); 

$uCase = length($block_cons);
#print "Ucase comes from $CONS_block\n";
##Not sure what to use here... if we go with CONS then we are including all the information found in the reads, even if they don't overlap with the sequences. 
$u_CONS = substr($CONS_block, $offset, $uCase); 
#print "$u_CONS\n\n";

$l_CONS = substr($CONS, ($uCase+1), length($CONS));
$l_CONS = lc $l_CONS;
#print "$l_CONS\n";

$CONS = $u_CONS.$l_CONS;

#print "Final cons is \n $CONS\n";


##NOW ADD THE SNPS to the CONS

@block = split('',$CONS);
#print Dumper (\@block);
foreach $snp (keys %hap_hash1){
	$pos = $snp;
	$pos = $pos -1; ##To adjust for zero start
	$pos = $pos - $offset; ###To adjust for multiple blocks
	#print "adding ". $hap_hash1{$snp}." to $pos\n";
	$block[$pos] = $hap_hash1{$snp};
}

$hap1 = join('',@block);
@block = split('',$CONS);
#print "$hap1 \n";


foreach $snp (keys %hap_hash2){
	$pos = $snp;
	$pos = $pos -1; ##To adjust for zero start
	$pos = $pos - $offset; ###To adjust for multiple blocks
	$block[$pos] = $hap_hash2{$snp};
}
$hap2 = join('',@block);

if($has_three == 1){
	@block = split('',$CONS);
	foreach $snp (keys %hap_hash3){
		$pos = $snp;
		$pos = $pos -1; ##To adjust for zero start
		$pos = $pos - $offset; ###To adjust for multiple blocks
		$block[$pos] = $hap_hash3{$snp};
	}
	$hap3 = join('',@block);
}

#$hap1=~s/N//g;
#$hap2=~s/N//g;
#$hap3 =~s/N//g;

print "OUTFILE IS: $outfile\n";
open FILE, ">", "$outfile";
print FILE ">$hap1_name\n$hap1\n>$hap2_name\n$hap2\n";
if($has_three==1){
	print FILE ">$hap3_name\n$hap3";
}
close FILE;

open (INPUT, $outfile) || die $!;
@edit_array = <INPUT>;
close(INPUT);

### REMOVE POSITONS THAT ARE ALL ? and - with ONLY ONE BASE! 
%gapPos = ();
$numSeq = 0;
foreach $line (@fas_array){
	if($line !~ />/){
		$numSeq = $numSeq + 1;
		$pos=0;
		#print "l comes from: ".$temp[0];
		$temp=$line;
		$temp = reverse($temp);
		#print "This is the seq $temp";
		while(length($temp)>0){
			$let=chop($temp);
			#if($let eq "-"){
			if($let eq "-" || $let eq "?"){
				$num = $gapPos{$pos};
				$num = $num +1;
				$gapPos{$pos} = $num;
			}
			$pos = $pos + 1;
		}
	}
}
$minus_one_numSeq = $numSeq - 1;

@deGapped = ();

print("DEGAPPING $outfile\n");
foreach $l (@edit_array){
	if($l !~/>/){
		$l = reverse($l);
		$currentPos=0;
		$newSeq ="";
		while(length($l)>0){
			$let=chop($l);
			#print "this is the key for $currentPos :".$seenPos{$currentPos};
			#if this position isn't a gap to be removed
			if($gapPos{$currentPos}==$numSeq || $gapPos{$currentPos} == $minus_one_numSeq){
			#print "$currentPos not included \n";
			}else{
				$newSeq=$newSeq.$let;
			}
			$currentPos = $currentPos + 1;
		}
		push(@deGapped, $newSeq);
	}
	else{
		push(@deGapped, $l);
	}
}
open FILE, ">" , "$outfile";
print FILE @deGapped;
close FILE;


###SUB ROUTINE####################

sub test3rd{
	@hapS1 = split('',$hapT1);
	@hapS2 = split('',$hapT2);
	@hapS3 = split('',$hapT3);
	
	
	$RECOMB = 0;
	##1 versus 2
	
	$l = scalar(@hapS1);
	
	$c = 1;
	
	
	#Comparing 1 and 2;
	while($c < $l && $RECOMB == 0){
		#switch each position
		#get the values at the position
		@pA = @hapS1[$c..$l];
		@pB = @hapS2[$c..$l];
		#print Dumper(@pA);
		#get a copy of the haplotypes
		@hapA = @hapS1;
		@hapB = @hapS2;
		$cm = $c - 1;
		#switch the positions
		@hapA=(@hapA[0..$cm],@pB);
		@hapB=(@hapB[0..$cm],@pA);
		
		#print "comparing ". join('',@hapA) ." and ".join('',@hapB)." with " . join('',@hapS3)."\n";
		
		if(join('',@hapA) eq join('',$hapS3) || join('',@hapB) eq join('',$hapS3)){
			$RECOMB = 3; 
		}
		$c = $c + 1;
	}
	$c = 1;
	@hapS1 = split('',$hapT1);
	@hapS2 = split('',$hapT2);
	@hapS3 = split('',$hapT3);
	
	#comparing 1 and 3
	
	while($c < $l && $RECOMB == 0){
		#switch each position
		#get the values at the position
		@pA = @hapS1[$c..$l];
		@pB = @hapS3[$c..$l];
		#print Dumper(@pA);
		#get a copy of the haplotypes
		@hapA = @hapS1;
		@hapB = @hapS3;
		$cm = $c - 1;
		#switch the positions
		@hapA=(@hapA[0..$cm],@pB);
		@hapB=(@hapB[0..$cm],@pA);
		
		#print "comparing ". join('',@hapA) ." and ".join('',@hapB)." with " . join('',@hapS2)."\n";
		
		if(join('',@hapA) eq join('',$hapS2) || join('',@hapB) eq join('',$hapS2)){
			$RECOMB = 2; 
		}
		$c = $c + 1;
	}
	
	$c = 1;
	@hapS1 = split('',$hapT1);
	@hapS2 = split('',$hapT2);
	@hapS3 = split('',$hapT3);
	
	#comparing 2 and 3
	
	while($c < $l && $RECOMB == 0){
		#switch each position
		#get the values at the position
		@pA = @hapS2[$c..$l];
		@pB = @hapS3[$c..$l];
		#print Dumper(@pA);
		#get a copy of the haplotypes
		@hapA = @hapS2;
		@hapB = @hapS3;
		$cm = $c - 1;
		#switch the positions
		@hapA=(@hapA[0..$cm],@pB);
		@hapB=(@hapB[0..$cm],@pA);
		
		#print "comparing ". join('',@hapA) ." and ".join('',@hapB)." with " . join('',@hapS1)."\n";
		
		if(join('',@hapA) eq join('',$hapS1) || join('',@hapB) eq join('',$hapS1)){
			$RECOMB = 1; 
		}
		$c = $c + 1;
	}
	
	
	
}
