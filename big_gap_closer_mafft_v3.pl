#!/usr/bin/perl
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

##V2 now accounts for MAFFT's innability to map fragments. It now does it in two steps: 
##requires teh split_frags.pl script

##V3 now doesn't overlap with 100 bp regions


#This file finds large gaps and extracts them to be analyzed with MUSCLE
##this step comes AFTER slider
if($#ARGV<0){
	print "*******************************************\nSyntax: big_gap_closer aligned_file.fasta";
	exit;
}



$file = $ARGV[0];

#open FILE, ">" , "$file.sites.all.fasta";
#close FILE;
	
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
		if($first == 1){
			$seq_len = length($line);
			#print "SEQ LEN IS $seq_len\n";
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

##Now we hav eot look for runs of high perc gaps 
$first = -1;
$first_found = -1;
$second = 0; 
$run = -1;

##NEED TO ACCOUNT FOR SQUISHING!! 
$adjustment = 0;

foreach $key (sort {$a<=>$b} keys %perc_gap){
	##CHECK ALL SPOTS!!!
	if($perc_gap{$key} > 0){
		#this might be the start or end
		if ($first == -1){
			#first time we have hit a gap
			$first = $key;
			$run = 0;
		}elsif($second == 0){
			$second = $key;
			$run = 0;
		}else{
			#STOP THE RUN IF WE HAVE GONE 100 BASES!!! then it turns into a sliding window. 
			if($run < 25 && ($second - $first) < 100){
				$second = $key;
				$run = 0;
				
			}else{
				#need to see if the next one is within 25? bases otherwise they will overlap in a difficult spot 
		
				#print "Checked upstream and found next_first to be $next_first and next_second to be $next_second\n"; 
				
				#print "RUN ENDED\n";
				$first = $first - $adjustment;
				$second = $second - $adjustment;
				print "Key is $key and First is : $first and Second is: $second and RUN is: $run \n";
				edit_muscle();
				#start over
				$first = -1;
				$second = 0;
				$run = 0;
			}
		}
	}else{
		#need to see if more than 8 have passed or else we need to RESET
		if($run == -1){
		}elsif($run > 25 && $second == 0){
			#start over
				$first = -1;
				$second = 0;
				$run = 0;
		}else{
			$run = $run + 1;
		}
	}
	
	
} 	

if($second != 0){
	
	#print "RUN ENDED at the end of the file\n";
	
	$first = $first - $adjustment;
	$second = $second - $adjustment;
	edit_muscle();  	
	print "Key is $key and First is : $first and Second is: $second and RUN is: $run \n";
}
	
sub edit_muscle{
	if(($second - $first) < 100){
		$start = $first - 25;
		if ($start <= 0){
			$start = 1;
		}
		$end = $second + 25;
	}else{
		$start = $first;
		$end = $second;
	}

	system("perl selectSites.pl -s '$start-$end' $file > sites.fasta");
	system("perl fasta_to_one_line.pl sites.fasta");
	open (INPUT, "sites.fasta");
	@original = <INPUT>;
	close(INPUT);
	system("perl split_frags.pl sites.fasta");
	if((-z "sites.fasta") == 1){
		#there were no complete sequences
		`mafft  --localpair --maxiterate 1000 --op 0 --lop 0 --nuc --quiet frags_sites.fasta > sites.out.fasta`;
	}else{
		`mafft  --localpair --maxiterate 1000 --op 0 --lop 0 --nuc --quiet sites.fasta > sites.out.fasta`;
	}
	
	if((-z "frags_sites.fasta") == 1){
		##there were no frags
	}else{
		`mafft --addfragments frags_sites.fasta --thread -1 --quiet sites.out.fasta > sites.fasta`;
		`cp sites.fasta sites.out.fasta`;
	}
		
	
	system("perl fasta_to_one_line.pl sites.out.fasta");
	#system("perl convert_terminal_gaps.pl sites.fasta");
	
	#now we have to add those back in! 
	open (INPUT, "sites.out.fasta");
	@fixed = <INPUT>;
	close(INPUT);
	
	
	
	
	
	#open FILE, ">>" , "$file.sites.all.fasta";
	#print FILE "\nKey is $key and First is : $first and Second is: $second and RUN is: $run \n";
	#print FILE "\nSelect sites is extracting from $start to $end\n";
	#print FILE @original;
	#close FILE;
	
	convert_terminal_gaps();
	
	#open FILE, ">>" , "$file.sites.all.fasta";
	#print FILE "\n\nFIXED\n";
	#print FILE @fixed;
	##close FILE;
	
	%fixed_hash = ();
	$current_id = "";
	foreach $line (@fixed){
		if($line =~ />/){
			$current_id = $line;
		}else{
			$line =~ s/\n//;
			$fixed_hash{$current_id} = $line;
		}
	}
	#adjust for 0 start
	$start = $start - 1;
	$end_pos = $end;
	$end = $end - $start;
	#SOMETIMES the blast output is not the same length as the input!!!!
	$fix_len = 0;
	foreach $f (@fixed){
		#print "$f\n";
		if($f !~ />/){
			$compare = length($f);
			#print "Comparing $compare to $end\n";
			if($compare != $end){
				$fix_len = 1;
				$to_add_adjust = $end - $compare;
				$adjustment = $adjustment + $to_add_adjust;
				$seq_len = $seq_len - $to_add_adjust;
				#print "adjustment = $adjustment\n and to_add_adjust = $to_add_adjust\n";
				#print "seq_len is now $seq_len\n";
				#print "LENTHS ARE DIFF \n";
			}
			last;
		}
		
	}
	
	
	$seq_len = $seq_len - $to_add_adjustment;
	#print "\nSelect sites is extracting from $start to $end_pos and adjusted len is: $seq_len\n";
	foreach $line (@fasta){
		if($line =~/>/){
			$current_id = $line;
		}else{
			if(defined $fixed_hash{$current_id} != 0){
				#print "Changing $line\n";
				$to_add = $fixed_hash{$current_id};
				#print "By adding $to_add\n";
				#print "to the positions $start and $end\n";
				##SUBSTR TAKES A START AND A LEN
				
				if($end_pos >= $seq_len){
					 substr($line, $start, $end) = ($fixed_hash{$current_id}."\n");
					 #print "NEED NEW LINE\n";
				}else{
					substr($line, $start, $end) = $fixed_hash{$current_id};
				}
				#print "To $line\n";
			}elsif($fix_len == 1){
				$to_add = "?"x$compare;
				substr($line, $start, $end) = $to_add;
			}
		}
	}
	
	#print Dumper (\%fixed_hash);
	
	#needs to be updated so that if i take positions out from the previous gap that it's alright
	open FILE, ">" , "$file";
	print FILE @fasta;
	close FILE;
#exit;
###FRONT PROBLEM!!! 0 START TURNS INTO -1 and BLAH BLAH
}
	
sub convert_terminal_gaps{
### This method needs to only change ones in @original that already had ??s at the start and finish. 	
	
##find the ones in @original that need to be changed. 
	%to_convert_front = ();
	%to_convert_end = ();
	$current_read = "";
	foreach $line (@original){
		#print "working with $line\n";
		if($line =~ />/){
			
			$current_read = $line; 
			#print "found a new read $current_read\n";
		}else{
			chomp($line);
			#print "This $current_read has $line as a line\n";
			$char = chop($line);
			#print "This $current_read has $char at the end\n";
			if($char eq "?"){
				$to_convert_end{$current_read} = 1;
				#print "This $current_read has ?s at the end\n";
			}
			
			$line = reverse($line);
			$char = chop($line);
			#print "This $current_read has $char at the front\n";
			if($char eq "?"){
				$to_convert_front{$current_read} = 1;
				#print "This $current_read has ?s at the front\n";
			} 
		}
	}
	$current_read = "";          
	foreach $line (@fixed){
		if($line =~ />/){
			$current_read = $line; 
		}
		@chars = ();
		if($line !~ />/){
			@chars = split('',$line);
			$l = scalar (@chars);
			$n = 0;
			if($to_convert_front{$current_read} == 1){				
				#print "This $current_read is being fixed at the front\n"; 
				while($n < $l){
					if($chars[$n] !~ /\-/){
						last
					}else{
						$chars[$n] = "?";
					}
					$n = $n + 1;
				}
				$n = $l;
				$n = $n-2;
			}
			$n = $l;
			$n = $n-2;
			if($to_convert_end{$current_read} == 1){
				#print "This $current_read is being fixed at the back\n"; 
				while($n > 0) {
					if($chars[$n] !~ /\-/){
						last;
					}else{
						$chars[$n] = "?";
					}
					$n = $n - 1; 
				}
			}
			$line = join('',@chars);
			$line =~ s/a/A/g;
			$line =~ s/t/T/g;
			$line =~ s/c/C/g;
			$line =~ s/g/G/g;
			
			
		}
	}
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
