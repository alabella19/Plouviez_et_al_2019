#!/usr/bin/perl
use Data::Dumper;
use Bio::AlignIO;
use Bio::SimpleAlign;
use List::Util qw(max);
use bignum;
$Data::Dumper::Sortkeys = 1;
####Finding the SNPS by looking down each column of a FASTA file

###V2. 8/6/13 can now handle FASTA files with sequences on multiple lines. 
####v3 8/7/13 will output the file needed to integrate with the R script.
####v6 10/20/14 Fixed the issue with SNPs not being called due to the lower snp perc being changed and NOT changed backs

#fasta file must be aligned. 
$file = $ARGV[0];
$upper_perc_snp = $ARGV[1];
$lower_perc_snp = $ARGV[2]; 
$indv = $ARGV[3];
$indv =~ s/\///;
$indv =~ s/\-/\./;	

$upper_perc_snp_saved = $upper_perc_snp;
$lower_perc_snp_saved = $lower_perc_snp; 

#print "PERC SNP = $upper_perc_snp";
if($upper_perc_snp == undef){
	$upper_perc_snp = 0.7;
}
if($lower_perc_snp == undef){
	$lower_perc_snp = 1 - $upper_perc_snp;
	$lower_perc_snp = $lower_perc_snp - .1 ; 
	#account for 10% error (there should be TWO and not 3 snp
}

#print "PERC SNP2 = $upper_perc_snp";
open(INPUT, $file) || die $!;
@input_array = <INPUT>;
close (INPUT);
@change_array = @input_array;

##Generate the consensus sequece from the original file
##in order to make a consensus sequence the ???s need to be turned back into ----
@new_array = ();
foreach $line (@change_array){
	if($line !~ />/){
		$line =~ s/\?/\-/g;
	}
	push(@new_array, $line);
}

open FILE, ">", "align_data.fasta";
print FILE @new_array;
close FILE;
	

$align_data = Bio::AlignIO->new(-file => "align_data.fasta");
$align_data = $align_data->next_aln();

$consensus_seq = $align_data->consensus_string(10);

#print "THIS IS CONSENSUS\n $consensus_seq\n";


@array_to_parse = @input_array;

##make a hash from each position in the fasta file 

%pos_hash = ();
%is_snp = ();
$seq = "";
$first = 0;
@fasta_array = ();
foreach $line (@array_to_parse){

	if($line !~ />/){
		$line =~ s/\n//;
		$seq = $seq.$line
	}elsif($first != 0){
		push(@fasta_array, $current_indv);
		push(@fasta_array, $seq);
		parse_line($seq);
		$seq = "";
		$current_indv = $line;
		$current_indv =~ s/\n//;
	}elsif($first == 0){
		$first = 1;
		$current_indv = $line;
		$current_indv =~ s/\n//;
	}
}
push(@fasta_array, $current_indv);
push(@fasta_array, $seq);
parse_line($seq);
$seq = "";
$current_indv = $line;
$current_indv =~ s/\n//;

###Now we need to assess if each position is a SNP or not??
##NUM SNPS
$num_snps = 0;
%is_N = ();
foreach $key (sort keys %pos_hash){
	##SNP = 0 means NOT A SNP
	$is_SNP_val = 0;
	$chars = $pos_hash{$key};
	
	@a = ($chars =~ /A/gi);
	$a = @a;
	@t = ($chars =~ /T/gi);
	$t = @t;
	@c = ($chars =~ /C/gi);
	$c = @c;
	@g = ($chars =~ /G/gi);
	$g = @g;
	@gap =($chars =~ /\-/gi); 
	$gap = @gap;
	$l = $a+$t+$c+$g+$gap;
	$readN = $l;
	if($l == 0){
		next;
	}
	
	##IF L is less than 10, then the minimum percentage for the minor allele should be more than 1 read"
	
	if($readN < 10){
		$lower_perc_snp = 1 - (($readN - 2)/$readN);
		$lower_perc_snp = sprintf("%.8f",$lower_perc_snp);
		#print "Lower SNP Percentage changed to : $lower_perc_snp for position $key\n";
	}else{
		#print "Lower SNP Percentage changed back to : $lower_perc_snp_saved for position $key\n";
		$lower_perc_snp = $lower_perc_snp_saved;
	}
	##HERE I HAVE COUNTED THE NUMBER OF READS THAT HAVE A SNP! IF THIS IS UNDER 4 WE WANT TO TURN IT INTO AN N
	#print "\tL is $l at $key";
		$a = $a/$l;
		$t = $t/$l;
		$c = $c/$l;
		$g = $g/$l;
		$gap = $gap/$l;
		$a = sprintf("%.8f",$a);
		$t = sprintf("%.8f",$t);
		$c = sprintf("%.8f",$c);
		$g = sprintf("%.8f",$g);
		$gap = sprintf("%.8f",$gap);
		
		
		@all_pos_list = ($a, $t, $c, $g, $gap);
		
		$max_perc = max @all_pos_list; 
		
		#check to see if anything is above the maximum SNP percentage. If it is then there is no SNP
		if($gap > $upper_perc_snp){
			$is_SNP_val = "gap"; 
		}elsif($a > $upper_perc_snp){
			$is_SNP_val = "A";
		}elsif($t > $upper_perc_snp){
			$is_SNP_val = "T";
		}elsif($c > $upper_perc_snp){
			$is_SNP_val = "C";
		}elsif($g > $upper_perc_snp){
			$is_SNP_val = "G";
		}else{
			#look for a minor alllele. Something greater than the lower perc snp but not = to the max
			#should this be > or >=
			if ($a != $max_perc && $a >= $lower_perc_snp){
				$is_SNP_val = 1;
				$num_snps = $num_snps + 1;
				#print "FOUND A MINOR ALLELE\n";
			}elsif($t != $max_perc && $t >= $lower_perc_snp){
				$is_SNP_val = 1;
				$num_snps = $num_snps + 1;
			}elsif($c != $max_perc && $c >= $lower_perc_snp){
				$is_SNP_val = 1;
				$num_snps = $num_snps + 1;
			}elsif($g != $max_perc && $g >= $lower_perc_snp){
				#print "ERROR OCCURS HERE  because  g = $g and lowerSNP = $lower_perc_snp\n";
				$is_SNP_val = 1;
				$num_snps = $num_snps + 1;
			}elsif($gap != $max_perc && $gap >= $lower_perc_snp){
				$is_SNP_val = 1;
				$num_snps = $num_snps + 1;
			}else{
				##there is no minor allele so go with the max perc base
				if($gap == $max_perc){
					$is_SNP_val = "gap"; 
				}elsif($a == $max_perc){
					$is_SNP_val = "A";
				}elsif($t == $max_perc){
					#something weird was going on here 11/5/14
					$is_SNP_val = "T";
				}elsif($c == $max_perc){
					$is_SNP_val = "C";
				}elsif($g == $max_perc){
					$is_SNP_val = "G";
				}
				
			}
			#need to account for the possibility of equal perecentages 
			if($a == $max_perc && $a > $lower_perc_snp){
				if($t == $max_perc || $c == $max_perc || $g == $max_perc || $gap == $max_perc){
					
					$is_SNP_val = 1;
					$num_snps = $num_snps + 1;
				}
			}elsif($t == $max_perc && $t > $lower_perc_snp){
				if($a == $max_perc || $c == $max_perc || $g == $max_perc || $gap == $max_perc){
					$is_SNP_val = 1;
					$num_snps = $num_snps + 1;
				}
			}elsif($c == $max_perc && $c > $lower_perc_snp){
				if($a == $max_perc || $t == $max_perc || $g == $max_perc || $gap == $max_perc){
					$is_SNP_val = 1;
					$num_snps = $num_snps + 1;
				}
			}elsif($g == $max_perc && $g > $lower_perc_snp){
				if($a == $max_perc || $c == $max_perc || $t == $max_perc || $gap == $max_perc){
					$is_SNP_val = 1;
					$num_snps = $num_snps + 1;
				}
			}elsif($gap == $max_perc && $gap > $lower_perc_snp){
				if($a == $max_perc || $c == $max_perc || $t == $max_perc || $g == $max_perc){
					$is_SNP_val = 1;
					$num_snps = $num_snps + 1;
				}
			}
			
			
			
			
			
			
		}
		
		if($readN <=5 && $is_SNP_val == 1){
			$is_SNP_val = 0; 
			$is_N{$key} = "N";
			#print "Is snp val is being set to 0 at $key\n";
			$num_snps = $num_snps-1;
		}else{
			$is_N{$key} = 0;
		}

		$is_snp{$key} = $is_SNP_val;
	
}
#print Dumper (\%is_N);
#print Dumper(\%pos_hash);
#print Dumper(\%is_snp);
#print Dumper(\%is_N);
#print $num_snps;
 	
$file_out = $file;

#print "FILE IS $file_out\n";
$file_out =~ s/\.fasta//;
$file_out =~ s/\.fna//;
$file_out =~ s/\.fas//;
$file_out_nex = $file_out.".haps.nex";
$file_consensus = $file_out.".consensus.txt";
$file_out = $file_out.".haps.fasta";


if($num_snps == 0){
	print "NO SNPS IN $file: CONSENSUS PRINTED to $file_out\n";
	print_simple();
}
if($num_snps == 1){
	print "ONLY 1 SNP IN $file: TWO HAPLOTYPES PRINTED to $file_out\n";
	print_simple();
}
if($num_snps >1){
	%nex_hash = ();
	$ntax = 0;
	print "MULTIPLE SNPS (N = $num_snps) IN $file :SNPS PRINTED TO $file_out_nex\n";
	@new_array = ();
	foreach $line (@fasta_array){
		if($line =~ />/){
			$current_read = $line;
			$current_read =~ s/>//;
			$current_read =~s/\n//;
			$current_read =~s/\r//;
			#print "cr = $current_read";
			$ntax = $ntax + 1;	
		}else{
			$new_line = get_snps($line);
			$nex_hash{$current_read} = $new_line;
		}
	}
	
	
	@nex_array = ();
	$taxa = "";
	$ntax=0;
	foreach $key (sort keys %nex_hash){
		$taxa = $taxa." ".$key;
		$ntax = $ntax+1;
	}
	
	push(@nex_array, "#NEXUS"."\t"x($num_snps)."\n[MacClade 4.08 registered to Department of Biology, Duke University]"."\t"x($num_snps-1)."\nBEGIN DATA;"."\t"x($num_snps-1));
	push(@nex_array, "\n\tDIMENSIONS NTAX=$ntax NCHAR=$num_snps;"."\t"x($num_snps-2));
	push(@nex_array, "\n\tFORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. INTERLEAVE;"."\t"x($num_snps-2));
	push(@nex_array, "\n\tCHARSTATELABELS"."\t"x($num_snps-2));
	$pos = "";
	#print Dumper \%is_snp;
	foreach $key (sort {$a<=>$b} keys %is_snp){
		#print "key is $key!\n";
		#print "value is $is_N{$key}\n";
		if($is_snp{$key} == 1 && $key > 0 && $is_N{$key} != 1){
			$pos = $pos.$key.", ";
			#print "key is $key!\n";
		}
	}
	
	push(@nex_array, "\n\t\t$pos"."\t"x($num_snps-3));
	push(@nex_array, "\nMATRIX"."\t"x($num_snps-1)."\n");
	foreach $key (sort {$a<=>$b} keys(%nex_hash)){
		#print "Key is $pos\n";
		$read = $nex_hash{$key};
		#print "$read\n";
		@read = split('',$read);
		$read = join('	',@read);
		#print "$read\n";
		$key =~ s/\n//;
		$key =~ s/\r//;
		$data = $key."	".$read."\n";
		push(@nex_array, $data);
	}
	push(@nex_array, ";\nEND;");
	open FILE, ">", "$file_out_nex";
	print FILE @nex_array;
	close FILE;
	
	#open FILE, ">", "$file_consensus";
	#print FILE "$pos\n$consensus_seq\n";
	#close FILE;
	
	##ALSO NEED TO GO BACK AND FOR EACH READ IN THE FILE ADD THE Ns in the RIGHT PLACES
	@with_N = ();
	#print $file;
	open(INPUT, $file) || die $!;
	@N_array = <INPUT>;
	close (INPUT);
	foreach $read (@N_array){
		if ($read =~ ">"){
			push(@with_N, $read);
		}else{
			@read = split('',$read);
			foreach $amb (keys %is_N){
				$N_val = $is_N{$amb};
				
				if($is_N[$amb] eq "N"){
					$amb = $amb -1; #to adjust for a zero start
					#print $read[$amb];
					print "isN value for key $amb is $N_val\n";
					if($read[$amb] !~ '\?'){
						$read[$amb] = "N";
					}
				}
			}
			$read = join('', @read);
			push(@with_N, $read);
		}
	}
	#print Dumper (\@with_N);
	open FILE, ">", "$file";
	print FILE @with_N;
	close FILE;
				
	
	
}

	
	


#####SUBROUTINES

sub get_snps{
	$seq = $_[0];
	$new_seq = "";
	foreach $pos (sort {$a<=>$b} keys %is_snp){
		if($is_snp{$pos} == 1){
			
			@seq_char = split("",$seq);
			$seq_char = $seq_char[($pos-1)];
			#print "Seqchar is $seq_char:";
			$seq_char=~ s/\?/N/;
			$new_seq = $new_seq.$seq_char;
		}
	}
	$new_seq;
}
	

sub parse_line{
	$seq = $_[0];
	$orig_l = length($seq);
	$seq =~ s/\n//g;
	$seq=~s/\r//g;
	while(length($seq) > 0){
			$l = length($seq);
			$char = chop($seq);
			#print "$l  ";
			$pos_chars = $pos_hash{$l};
			$new_pos_chars = $pos_chars.$char;
			$pos_hash{$l} = $new_pos_chars;
	}
}

sub print_simple{
	##To add ambiguities and then add the 1 SNP if there is one.
	#what is the threshold for assigning an ambiguity 
	#IUPAC doesn't account for - versus ATGC 
	#consensus_string counts ?s in the %
	open(INPUT, $file) || die $!;
	@q_array = <INPUT>;
	close (INPUT);

	%rev_hash = reverse %is_snp;
	@s_snps = ();
	$single_snp = -1;
	#print keys %rev_hash;
	if(exists $rev_hash{1}){
		
		$single_snp = $rev_hash{1};
		print "THERE IS A SNP at $single_snp\n";
		#print "pos hash?";
		#print $pos_hash{$single_snp};
		
		$chars = $pos_hash{$single_snp};
		
		@a = ($chars =~ /A/gi);
		$a = @a;
		$chars{"A"} = $a;
		@t = ($chars =~ /T/gi);
		$t = @t;
		$chars{"T"} = $t;
		@c = ($chars =~ /C/gi);
		$c = @c;
		$chars{"C"} = $c;
		@g = ($chars =~ /G/gi);
		$g = @g;
		$chars{"G"} = $g;
		@gap =($chars =~ /\-/gi); 
		$gap = @gap;
		$chars{"-"} = $gap;
		 %chars = reverse %chars;
		#print Dumper (\%chars);
		##Just get the top two bases
		
		#print max keys %chars;
		$base_1 = $chars{(max keys %chars)};
		$base_1_val = max keys %chars;
		delete $chars{(max keys %chars)};
		#print Dumper (\%chars);
		$base_2 = $chars{(max keys %chars)};
		$base_2_val = max keys %chars;
		#print "base 1 is $base_1\n";
		#print "base 2 is $base_2\n";
		
		###NEED TO ADD IN THE BASES NOW!
	}

	foreach $line (@q_array){
			$line =~ s/-/N/g;
			#$line =~ s/\?/-/g;
		}
	open FILE, ">", "cons.fasta";
	print FILE @q_array;
	close FILE;
	
	$in = Bio::AlignIO->new(-file => 'cons.fasta' , -format => 'fasta');
	$align = $in->next_aln();
	$cons_iupac_N = $align->consensus_iupac();
	
	
	foreach $line (@q_array){
			$line =~ s/N/-/g;
			$line =~ s/\?/-/g;
	}
	
	open FILE, ">", "cons.fasta";
	print FILE @q_array;
	close FILE;
	
	$in = Bio::AlignIO->new(-file => 'cons.fasta' , -format => 'fasta');
	$align = $in->next_aln();
	$cons_iupac = $align->consensus_iupac();
	$cons_thresh = $align->consensus_string(30);
	
	#print Dumper(\%is_snp);
	#print Dumper(\@thresh);
	#print Dumper(\@iupac);
	#value will be 1 for a key if it is a SNP
	##get IUPAC N
		
	
	
	#print("IUPAC Consensus is \n$cons_iupac\n and thresh is \n$cons_thresh\n and iupac_N is \n$cons_iupac_N");
	@iupac = split('',$cons_iupac);
	@thresh = split ('',$cons_thresh);
	@iupac_N = split ('',$cons_iupac_N);
	@final_cons = ();
	
	
	$pos = 0; 
	foreach $i_base (@iupac){
		#print("comparing $thresh[$pos] and $iupac[$pos]\n");
		if($thresh[$pos] =~ /\?/ && $iupac[$pos] !~ /[AaTtCcGg]/ && $iupac_N[$pos] !~ /n/){
			#if nothing meets the threshold go with iupac unless iupac_N is n (aka there is a gap)
			push(@final_cons, uc($iupac[$pos]));
		}elsif($thresh[$pos] =~ /\?/ && $iupac[$pos] =~ /[AaTtCcGg]/ && $iupac_N[$pos] !~ /n/){
			#if nothing meets the threshold then just go with iupac unless iupac_N is n (aka there is a gap)
			push(@final_cons, uc($iupac[$pos]));
		}elsif($thresh[$pos] !~ /\?/ && $iupac[$pos] !~ /[AaTtCcGg]/){
			#if something meets the threshold and there is an ambiguity, go with the thresh
			push(@final_cons, uc($thresh[$pos]));
		}elsif($thresh[$pos] !~ /\?/ && $iupac[$pos] =~ /[AaTtCcGg]/){
			#if somsething meets the threshold, go with that;
			push(@final_cons, uc($thresh[$pos]));
		}elsif($thresh[$pos] =~ /\?/ && $iupac_N[$pos] =~ /n/){
			#if there is a gap go with LOWERCASE iupac
			#print "found a possible gap at $pos\n";
			push(@final_cons, $iupac[$pos]);
		}
		

		#print"pos is $pos\n";
		$pos = $pos + 1;
		
	}

	$final_cons = join('',@final_cons);
	#print "Final cons is $final_cons\n\n";
	
	
	#print Dumper (\%is_snp);
	#add in the gaps 
	foreach $pos (keys %is_snp){
		if($is_snp{$pos} =~ 'gap'){
			substr($final_cons, ($pos - 1), 1) = "-";
		}
	}
	
	#now we need to add in $base_1 and $base_2 to $single_snp
	if($single_snp != -1){
		$half_1 = substr($final_cons, 0, $single_snp);
		$half_2 = substr($final_cons, $single_snp, length($final_cons));
		
		##SNP is on half 1
		chop($half_1);
		$half_1a = $half_1 . $base_1;
		$half_1b = $half_1 . $base_2;
		$hap1=$half_1a . $half_2;
		$hap2=$half_1b . $half_2;
		#print "$hap1\n\n";
		#print "$hap2\n"; 
	}else{
		$hap1 = $final_cons;
		$hap2 = $final_cons;
		$base_1_val = "hom";
		$base_2_val = "hom";
	}
	
	open FILE, ">", "T_$file_out";
	print FILE ">$indv.1.$base_1_val\n$hap1\n>$indv.2.$base_2_val\n$hap2";
	close FILE;
}

sub snp_val() {
    $large_key = -1;
    $large_val = -1;
    foreach $key (keys %val_hash) {
    	    $val = $val_hash{$key};
    	    #print "Key is $key\n";
        if ($val > $large_val) {
            $large_val = $val;
            $large_key = $key;
        }
    }
    $large_key
}
