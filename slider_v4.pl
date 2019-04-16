#!/usr/bin/perl
use Data::Dumper;


##version 3 updated on 10.24.14 to account for -- in the squishing.
##version 4 checks to make sure the large portion is the majority of the reads. Otherwise it doesn't reduce the size
#This file collapses or slides together regions that are screwed up because of a low frequency indel
##Example:
##	A T - G C G T
##	A T - G C G T
##	A T G C G - T
##	A T G C G - T
##      A T G C G G T
##
##	Would turn into 

##	A T G C G T
##	A T G C G T
##	A T G C G T
##	A T G C G T
##      A T G C G T
##
##I really don't want to mess with anything 10 BP or longer... 
## Or with more than 3 gaps
if($#ARGV<0){
	print "*******************************************\nSyntax: slider.pl aligned_file.fasta";
	exit;
}

$file = $ARGV[0];

open (INPUT, $file);
@fasta = <INPUT>;
close(INPUT);

%gapPos = ();
$numSeq = 0;
foreach $line (@fasta){
	if($line !~ />/){
		$numSeq = $numSeq + 1;
		$pos=0;
		#print "l comes from: ".$temp[0];
		$temp=$line;
		$temp = reverse($temp);
		#print "This is the seq $temp";
		while(length($temp)>0){
			$let=chop($temp);
			if($let eq "-"){
			#if($let eq "=" || $let eq "N"){
				$num = $gapPos{$pos};
				$num = $num +1;
				$gapPos{$pos} = $num;
			}
			$pos = $pos + 1;
		}
	}
}

#print Dumper(\%gapPos);	
#print "Total number of seqs on $outfile is $numSeq\n";
@deGapped = ();

print("DEGAPPING $outfile\n");
foreach $l (@fasta){
	if($l !~/>/){
		$l = reverse($l);
		$currentPos=0;
		$newSeq ="";
		while(length($l)>0){
			$let=chop($l);
			#print "this is the key for $currentPos :".$seenPos{$currentPos};
			#if this position isn't a gap to be removed
			if($gapPos{$currentPos}==$numSeq){
			#print "$currentPos not included ";
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

#REMOVE ALL OF THE POSITONS THAT ONLY HAVE ONE BASE
%gapPos = ();
$numSeq = 0;
foreach $line (@deGapped){
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




#print Dumper(\%gapPos);	
#print "Total number of seqs on $outfile is $numSeq\n";

#allows for a single base in the column -> column will be removed
$minus_one_numSeq = $numSeq - 1;

@deIndel = ();



print("DEGAPPING $outfile\n");
foreach $l (@deGapped){
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
		push(@deIndel, $newSeq);
	}
	else{
		push(@deIndel, $l);
	}
}

@fasta = @deIndel;
 



#get the positions of the rows that have 30% Gap or more 
%to_be_fixed = ();
$tot_to_be_fixed = 0;
%gapPos = ();
%tot_reads = ();
$numSeq = 0;
$first = 1; 
foreach $line (@fasta){
	if($line !~ />/){
		$numSeq = $numSeq + 1;
		$pos=0;
		#print "l comes from: ".$temp[0];
		$temp=$line;
		$temp = reverse($temp);
		#print "This is the seq $temp";
		while(length($temp)>0){
			
			if($first == 1){
				$gapPos{$pos} = 0; 
			}
			
			$let=chop($temp);

			if($let eq "-"){
				$num = $gapPos{$pos};
				$num = $num +1;
				$gapPos{$pos} = $num;
			}
			if($let ne "?"){
				$num_r = $tot_reads{$pos};
				$num_r = $num_r +1;
				$tot_reads{$pos} = $num_r;
			}
			$pos = $pos + 1;
		}
	$first =0;
	}
}
$tot_pos = $pos - 1; 
%perc_gap = ();


foreach $key (keys %gapPos){
	#sprint $key;
	$gaps = $gapPos{$key};
	$tot = $tot_reads{$key};
	$perc = $gaps / $tot;
	$perc_gap{$key} = $perc;
}

#print Dumper (\%perc_gap);

#exit;
$first = -1;
$first_found = -1;
$second = 0; 
foreach $key (sort {$a<=>$b} keys %perc_gap){
	#.1 might have been too liberal. .2 might be better
	if($perc_gap{$key} >.2 && $perc_gap{$key} <.8){
		#print "Found gap $key\n";
		if ($first == -1){
			$first = $key;
		}elsif($first != -1 && $second == 0){
			$second = $key;
			#print "Bounds are $first and $second\n";
			$result = compare_reads();
			if($result == 1){
				#worked and need to move on
				$first = -1;
				$second = 0;
			}
		}else{
			$first = $second;
			$second = $key;
			#print "Bounds are $first and $second\n";
			$result = compare_reads();
			if($result == 1){
				#worked and need to move on
				$first = -1;
				$second = 0;
			}
		}

		
	}
	
}

###SUBROUTINE

sub compare_reads{
	$read_num = 0;
	%reads = ();
	$lower = $first; # - 1;
	$length = $second - $lower + 1 ;
	$max_len = 0;
	$min_len = -1;
	%lengths = (); 
	$num_diff_len = 0;
	$first_zero = 0;
	#print Dumper (\@fasta);
	#$q_mark = "?" x $length;
	#print "Trying to make a substr starting at $lower that is $length long and compared to $q_mark\n";
	#if there are more than 2 lengths then this won't work and it should be aborted
	foreach $line (@fasta){
		if($line !~ />/){
			$this_read = substr($line, $lower, $length);
			if ($this_read !~ /\?/){
				#print "this_read is $this_read\n";
				$this_read =~ s/-//g;
				$reads{$read_num} = $this_read;
				
				if($this_read eq ""){
					#print "THERE IS A ZERO LENGTH\n";
					$this_len = "zero";
				}else{
					$this_len = length($this_read);
				}
				##I only want this to work if there are 2 different lengths, not 3 since that gets super confusing
				if(defined $lengths{$this_len} == 0){
					if($this_len eq "zero"&& $first_zero == 0){
						#print "ignoring first zero\n";
						
						#if this is the first time we hav ehad a zero length then ignore
					}else{
						$num_diff_len = $num_diff_len + 1;
						$lengths{$this_len} = 1; 
					}
					#print "FOUND A NEW LENGTH\n";
				}
				#print"NUM DIFF LENGTHS\n";
				#print Dumper (\%lengths);
				if($num_diff_len == 3){
					#print "Three Different Lengths and RETURN\n";
					return 0;
				}
				
				if($this_len > $max_len){
					$max_len = $this_len;
				}
				#print "This len is $this_len\n";
				if($this_len < $min_len || $min_len == -1){
					if($this_len ne "zero"){
						#print "setting min len to $this_len\n";
						$min_len = $this_len;
					}
				}
				if($this_len eq "zero"){
					#print "Min set to zero\n";
					if($first_zero == 0){
						$first_zero = 1;
					}else{
						#print "found second zero!\n";
						$min_len = 0;
					}
				}
				
				
				$read_num = $read_num + 1;
			}
		}
	}
	$len_diff = $max_len - $min_len;
	$num_gaps = $length - $max_len;
	#print "Len diff $len_diff = Max Len $max_len minus Min Len $min_len and total len is: $length\n";
	if($len_diff ==1 && $length < 11 && $num_gaps < 3){
		remove_split();
		#WORKED and these positions shouldn't be used again
		###COULD ADD LENGHT DIFFERENCE = 0 THEN YOU WOULD REMOVE THE SPACES 
		return 1;
	}elsif($len_diff == 0 && $length < 11){
		$gaps_removed = remove_gaps();
		if($gaps_removed == 1){
			return 1;
		}else{
			return 0;
		}
	}else{
		#DIDN'T WORK
		return 0;
	}
	
	#print Dumper (\%reads);
	#print "Min length is $min_len and max length is $max_len";
}
sub remove_gaps{
	#returns zero if they weren't all the same and 1 if they were
	#In order for this to work they must all be the same
	print "REMOVING GAP_Split\n";
	
	#make sure they are all teh same first
	#print "ALL THESE MUST BE THE SAME\n";
	#print Dumper (\%reads);
	@keys = keys %reads; 
	$key_first = $keys[0];
	$match = 0;
	foreach $key_check (@keys){
		if($reads{$key_check} ne $reads{$key_first}){
			#print "THEY DID NOT ALL MATCH!\n";
			##Could this be a SNP? Check if it is one bp long?
			$match = 1; 
			last;
		}
	}
	if($match == 1){
		#print "Is it a SNP?\n";
		if($min_len == $max_len){
			#Check to see if it is a SNP
			%SNP = ();
			foreach $read (keys %reads){
				if(defined $SNP{$reads{$read}} == 0){
					$SNP{$reads{$read}} = 1;
				}
			}
			#print "I THINK THERE MIGHT BE A SNP\n";
			#print Dumper (\%SNP);
			$len_SNP = scalar (keys %SNP);
			#print "NUM SNPS IS $len_SNP\n";
			if($len_SNP != 2){
				return 0;
			}
			
		}else{
			return 0;
			#print "Not a SNP :-(\n";
		}
	}
	#if they all match or the length is 1 and there are 2 bases you can proceed
	print "SQUISHING \n";
	$to_add_length = $length - $min_len;
	$to_add = "-" x $to_add_length;
	foreach $line (@fasta){
		
		if($line !~ />/){
			$this_read = substr($line, $lower, $length);
			if ($this_read !~ /\?/){
				#print "Replacing $this_read ";
				$this_read =~ s/-//g;
				if($this_read eq ""){
					$this_read = "-" x $length;
				}else{
					$this_read = $this_read.$to_add;
				}
				#print "With this: $this_read \n";
				substr($line, $lower, $length) = $this_read;
			}
			else{
				$replace = "?" x $length;
				substr($line, $lower, $length) = $this_read;
			}
		}
	}
	return 1; 
}

sub remove_split{
	print "REMOVING SPLIT\n";
	#gap gets added to the end of each. 
	$to_add_length = $length - $min_len;
	$to_add = "-" x $to_add_length;
	$to_add_length_max = $length - $max_len;
	$to_add_max = "-" x $to_add_length_max;
	%match_min = ();
	
	##NEED TO SEE IF IT IS A LINNKED GAP that was misaligned
	%match_all = (); 
	%match_all_num = (); 
	foreach $line (@fasta){
		if($line !~ />/){
			$this_read = substr($line, $lower, $length);
			if ($this_read !~ /\?/){
				$this_read =~ s/-//g;
				if($max_len != $min_len && length($this_read) == $max_len){
					$this_read = $this_read.$to_add_max;
				}
				$match_num = $match_all{$this_read};
				$match_num = $match_num + 1;
				$match_all{$this_read} = $match_num; 
				$match_all_num{$this_read} = $match_num;
				
			}
		}
	}
	
	#print "Need to check on linked gap\n";
	#print Dumper (\%match_all);
	if(scalar keys %match_all == 2){
		##Could be a linked SNP
		$total = 0;
		foreach $key (keys %match_all){
			$total = $total + $match_all{$key};
		}
		foreach $key (keys %match_all){
			$num = $match_all{$key};
			$num = $num/$total;
			if($num < .6 && $num > .4){
				#print "PROBABLY A LINKED GAP!";
				##and we want to put the gaps all in teh same place
				foreach $line (@fasta){
					if($line !~ />/){
						$this_read = substr($line, $lower, $length);
						if ($this_read !~ /\?/){
							$this_read =~ s/-//g;
							$reads{$read_num} = $this_read;
							$this_len = length($this_read);
							if($this_len == $min_len){
								$match_min{$this_read} = 1;
								$this_read = $this_read.$to_add;
								substr($line, $lower, $length) = $this_read;
							}
						}
					}
				}
			return;
			}
		}
	
	}
	
	##So now I need to see where to add teh gap! check in match all
	if($min_len > 1){
		#check to see what to put for the little ones. 
		$to_compare = "";
		$found = 0;
		foreach $key (keys %match_all){
			if(length ($key) == $min_len && $found == 0){
				#print "Setting compare\n";
					#if it is a short one and there isn't one to compare
					$to_compare = $key . $to_add;
					@to_compare_char = split(//,$to_compare);
					$match_all{$key} = $to_compare;
					$found = 1;
			}elsif(length ($key) == $min_len){
				$front_score = 0;
				$back_score = 0;
				##figure out which of the two is better;
				$test_front = $to_add . $key;
				$test_back = $key . $to_add; 
				@test_front_char = split(//,$test_front);
				@test_back_char = split(//,$test_back);
				
				#print "Testing $test_front versus $test_back \n";
				$p = 0;
				while($p<length @to_compare_char){
				#	print "Comparing front: " . $test_front_char[$p] . " and back: " . $test_back_char[$p];
				#	print " With " . $to_compare_char[$p] . "\n";
					if($test_front_char[$p] ne $to_compare_char[$p]){
						$front_score = $front_score + 1;
					}
					if($test_back_char[$p] ne $to_compare_char[$p]){
						$back_score = $front_score + 1;
					}
					$p = $p +1;
				}
				#print "scores are front: $front_score and back: $back_score\n";
				if($front_score < $back_score){
					#print "FRONT is better: $test_front\n";
					$match_all{$key} = $test_front;
				}elsif($back_score < $front_score || $front_score == $back_score){
					$match_all{$key} = $test_back;
					#print "Back is better $test_back\n"
				}
				
			}
		}
	}
	
	#print "New match all has is:\n";
	#print Dumper (\%match_all);
	

	
	foreach $line (@fasta){
		if($line !~ />/){
			$this_read = substr($line, $lower, $length);
			if ($this_read !~ /\?/){
				#DO I WANT THIS: DOES THIS TRHOW OFF THE READS!!!
				$this_read =~ s/-//g;
				$reads{$read_num} = $this_read;
				$this_len = length($this_read);
				if($this_len == $min_len){
					$match_min{$this_read} = 1;
					if($min_len > 1){
						$this_read = $match_all{$this_read};
					}else{
						$this_read = $this_read.$to_add;
					}
					substr($line, $lower, $length) = $this_read;
				}elsif($this_len == $max_len && $max_len != $len){
					$this_read = $this_read.$to_add_max;
					substr($line, $lower, $length) = $this_read;
				}

				#$read_num = $read_num + 1;
			}else{
				$replace = "?" x $length;
				substr($line, $lower, $length) = $this_read;
			}
		}
	}
	
	
	#Only go onto the next step if the gap is the prominent read
	#print Dumper (\%match_all_num);
	#print "to match are ";
	#print Dumper (\%match_min);
	#print "\n";
	$sum_min;
	$sum_max;
	#print "min len is $min_len and max len is $max_len\n";
	foreach $key (keys %match_all_num){
		
		if(length($key) >= $max_len){
			#print "key is max $key and length is" . length($key);
			$sum_max = $sum_max + $match_all_num{$key};
		}elsif(length($key) == $min_len){
			#print "key is min $key and length is" . length($key);
			$sum_min = $sum_min + $match_all_num{$key};
		}
	}
	
	$min_perc_reads = $sum_min/($sum_min + $sum_max);
	
	if($min_perc_reads>0.9){
		print("Removing positions because min reads make up $min_perc_reads of the file\n");
		foreach $line (@fasta){
			if($line !~ />/){
				$this_read = substr($line, $lower, $length);
				if ($this_read !~ /\?/){
					$this_read =~ s/-//g;
					$reads{$read_num} = $this_read;
					$this_len = length($this_read);
					if($this_len == $max_len){
						%matched_bases = ();
						#print "found a long one $this_read\n";
						#which base in the one that is too long is incorrect 
						$num_tries = $length - $min_len;
						$offset = 0;
						while($offset <= $num_tries){
							$test = substr($this_read, $offset, $min_len);
							#print "$test \n";
							#which test is the correct test for this read
							if (defined $match_min{$test} != 0){
								#print "found a match!\n";
								$matched_bases{$test} = $test;
							}
							$offset = $offset + 1;
						}
						#print "Matched bases are";
						#print Dumper (\%matched_bases);
						#print "\n";
						if(scalar keys %matched_bases == 1){
							#There is only one subset that matches a known allele
							@found = keys %matched_bases; 
							$found = $found[0];
							$found = $found.$to_add;
							substr($line, $lower, $length) = $found;
							print "\nFound a replacement for the insertion: $this_read replaced with $found\n";
						}else{
							##There is more than one match and we cannot determine what this should be
							print "\nCould not find a replacement for the insertion: $this_read\n";
							$test = "N" x $min_len;
							$test_to_add = "-" x $to_add_length;
							$test = $test.$test_to_add;
							substr($line, $lower, $length) = $test;
						}
					}
	
					#$read_num = $read_num + 1;
				}
			}
		}
	}
}



#Now remove all of the posit ions that are all gaps!

#get rid of the positions that are ALL gaps
%gapPos = ();
$numSeq = 0;
foreach $line (@fasta){
	if($line !~ />/){
		$numSeq = $numSeq + 1;
		$pos=0;
		##print "l comes from: ".$temp[0];
		$temp=$line;
		$temp = reverse($temp);
		#print "This is the seq $temp";
		while(length($temp)>0){
			$let=chop($temp);
			if($let eq "-"){
			#if($let eq "=" || $let eq "N"){
				$num = $gapPos{$pos};
				$num = $num +1;
				$gapPos{$pos} = $num;
			}
			$pos = $pos + 1;
		}
	}
}

#print Dumper(\%gapPos);	
#print "Total number of seqs on $outfile is $numSeq\n";
@deGapped = ();

print("DEGAPPING $outfile\n");
foreach $l (@fasta){
	if($l !~/>/){
		$l = reverse($l);
		$currentPos=0;
		$newSeq ="";
		while(length($l)>0){
			$let=chop($l);
			#print "this is the key for $currentPos :".$seenPos{$currentPos};
			#if this position isn't a gap to be removed
			if($gapPos{$currentPos}==$numSeq){
			#print "$currentPos not included ";
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

#REMOVE ALL OF THE POSITONS THAT ONLY HAVE ONE BASE
%gapPos = ();
$numSeq = 0;
foreach $line (@deGapped){
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

#print Dumper(\%gapPos);	
#print "Total number of seqs on $outfile is $numSeq\n";

#allows for a single base in the column -> column will be removed
$minus_one_numSeq = $numSeq - 1;

@deIndel = ();

print("DEGAPPING $outfile\n");
foreach $l (@deGapped){
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
		push(@deIndel, $newSeq);
	}
	else{
		push(@deIndel, $l);
	}
}



open FILE, ">" , "$file";
print FILE @deIndel;
close FILE;
#print Dumper (\%perc_gap);

