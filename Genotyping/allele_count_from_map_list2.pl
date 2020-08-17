use strict;
use warnings;
use Class::Struct;
use Getopt::Std;
use IO::File;
use Time::HiRes;
use Cwd;
my % options=();
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Handle;

my$opt_string='hm:o:p:g:n:r:';

getopts("$opt_string",\%options) or usage();
usage()if $options{h};
#change here if you have a gzip file
my $GZIP=0;

my $REPORT_ALL = 0;


my $out_file_handler = IO::File->new();
my $mapDotList="";
my %POS = ();
my %Cons=();
my %POS_ref=();
my $output_file_name="output.txt";
my $support=0;
my %length_distr=();
my $marker=0;
my $name=0;
my $pattern="pattern.txt";
my $max=100000000000;
my $restrictions_site_file="";
my %restrictions_site_rs_cons=();
my %restrictions_site_rs_cons_map=();
my $cut_length=5;
my $new_map_list_file="new_filtered_map.list";
my %quality_base_value=();

init();
readMarkers();
read_rs_sites();
main();

sub usage
{
	print "perl $0 -n name -m map.list -p marker -r restrictions_sites -o output -g beginning and end cut\n";
}

sub init
{
	if($options{m} && $options{p} && $options{n} && $options{r})
	{
		$mapDotList=$options{m};
		$marker=$options{p};
		$name=$options{n};
		$restrictions_site_file = $options{r};
		
	}
	else
	{
		die usage();
	}	
	if($options{o})
	{
		$output_file_name=$options{o};
	}
	if($options{g})
	{
		$cut_length=$options{g};
	}
}



sub read_rs_sites
{
	my $input_file = new IO::File($restrictions_site_file, "r") or die "could not open $restrictions_site_file: $!\n";
	while (my $line = $input_file->getline) 
	{
		chomp($line);
		my @a = split " ",$line;
		
		if(!exists $restrictions_site_rs_cons{$a[0]." ".$a[1]})
		{
			my @temp = ("A","A");
			$restrictions_site_rs_cons_map{$a[0]." ".$a[1]}=[@temp];
#			if($a[0] == 1 && $a[1]==371052)
#			{
#				my @bla = @{$restrictions_site_rs_cons_map{$a[0]." ".$a[1]}}; 
#				print "@bla \n";
#			}
			$restrictions_site_rs_cons{$a[0]." ".($a[1]-10)}=$a[0]." ".$a[1];
			$restrictions_site_rs_cons{$a[0]." ".$a[1]}=$a[0]." ".$a[1];
			$restrictions_site_rs_cons{$a[0]." ".($a[1]+10)}=$a[0]." ".$a[1];		
		}
	}
}


sub update_snp_counter
{ 
	my ($p1_read_array_ref,$p2_read_array_ref,$hash_value)=@_;
	my ($p1_beg,$p1_end) = update_snp_under_function($p1_read_array_ref); 
	my ($p2_beg,$p2_end) = update_snp_under_function($p2_read_array_ref);
	
	my @a=();
	#print " hash_value ".$hash_value."\n";
	my @temp = @{$restrictions_site_rs_cons_map{$hash_value}};
	#print " @temp\n";
	if(@{$restrictions_site_rs_cons_map{$hash_value}}[0] eq "A")
	{
		if($p1_beg  != -1)
		{
			push(@a,$p1_beg);
		}
		if($p1_end != -1)
		{
			push(@a,$p1_end);
		}
		if($p2_beg  != -1)
		{
			push(@a,$p2_beg);
		}
		if($p2_end !=-1)
		{
			push(@a,$p2_end);
		}	
	}
	else
	{
		if($p1_beg  != -1)
		{
			push(@a,$p1_beg);
		}
		if($p1_end != -1)
		{
			push(@a,$p1_end);
		}
		if($p2_beg  != -1)
		{
			push(@a,$p2_beg);
		}
		if($p2_end !=-1)
		{
			push(@a,$p2_end);
		}		
		push(@a,$restrictions_site_rs_cons_map{$hash_value}[0]);
		push(@a,$restrictions_site_rs_cons_map{$hash_value}[1]);	
	}
	#now sort it !
	if(@a >0)
	{
		my @sorted_values = sort {$a<=>$b} @a;
		$restrictions_site_rs_cons_map{$hash_value}=[$sorted_values[0],$sorted_values[@sorted_values-1]];
	}
}


sub update_snp_under_function
{
	my ($read_array_ref)=@_;
	my @a = @$read_array_ref;
	
	#print " array @a \n";	
	my $chr=$a[0];
	my $read_seq_length=$a[7];
	my $sequence=$a[2];
	my $pos_counter=$a[1];
	my $begin_ = -1;
	my $end_ = -1;
	my $i=0;
	my %tempMarkerFile=();
	my @tempMarkerVal=();
	my $read_length_counter=0;
	my $quality_string = $a[10];
	if($a[4] eq "P")
	{
		$quality_string = scalar reverse($quality_string); 
	}
#go through the sequence;
#	print($chr." ".$pos_counter."\n");
#	print($sequence."\n");
	while($i<length($sequence))
	{
		#print(substr($sequence,$i,1)." ");
		#print" pos: ".$pos_counter."\n";
		if(substr($sequence,$i,1) eq "[") #indel are not considered as a marker for our case!
		{
		#	print"mismatch or deletion\n";
		#	print" pos: ".$pos_counter."\n";
			++$i;
		#	print(substr($sequence,$i,1)." ");
		#	print" pos: ".$pos_counter."\n";
			my $baseREF=substr($sequence,$i,1);
			++$i;
		#	print(substr($sequence,$i,1)." ");
		#	print" pos: ".$pos_counter."\n";
			my $baseSNP=substr($sequence,$i,1);
			++$i;
		#	print(substr($sequence,$i,1)." ");
		#	print" pos: ".$pos_counter."\n";
			if($baseREF ne "-" && $baseSNP ne "-")
			{
				#$a[1]."#".$a[2];
				if(exists $POS{$chr."#".$pos_counter} && $read_length_counter>=$cut_length && $read_length_counter<$read_seq_length-$cut_length )
				{
					if($begin_ == -1)
					{
						$begin_ = $pos_counter;
					}
					else
					{
						$end_ = $pos_counter;
					}
					addQualityValue($quality_string,$read_length_counter);
					$tempMarkerFile{$chr."#".$pos_counter}=[0,0,0,0];
					my $name=$chr."#".$pos_counter;
		#			print"init"."@{$tempMarkerFile{$name}}"."\n";
					push(@tempMarkerVal,$chr."#".$pos_counter);
					
					#update counts
					if($baseSNP eq "A")
					{
						my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
						$temp[0]+=1;
						$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
		#				print"after a "."@{$tempMarkerFile{$name}}"."\n";		
					}
					elsif($baseSNP eq "C")
					{
						my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
						$temp[1]+=1;
						$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
		#				print"after c "."@{$tempMarkerFile{$name}}"."\n";
					}
					elsif($baseSNP eq "G")
					{
						my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
						$temp[2]+=1;
						$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
		#				print"after g "."@{$tempMarkerFile{$name}}"."\n";
					}
					elsif($baseSNP eq "T")
					{
						my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
						$temp[3]+=1;
						$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
		#				print"after t "."@{$tempMarkerFile{$name}}"."\n";
					}
		#			print"after"."@{$tempMarkerFile{$name}}"."\n";							
				}
				$pos_counter++;
				++$read_length_counter;
				
			}
			else 
			{
				if ($baseREF ne "-") 
				{
					$pos_counter++;
				}
				else
				{
					++$read_length_counter;
				} 
			}
		}
		else
		{
			if(exists $POS{$chr."#".$pos_counter} && $read_length_counter>$cut_length && $read_length_counter<$read_seq_length-$cut_length)
			{
				if($begin_ == -1)
				{
					$begin_ = $pos_counter;
				}
				else
				{
					$end_ = $pos_counter;
				}
				addQualityValue($quality_string,$read_length_counter);
				my $baseREF=substr($sequence,$i,1);
				$tempMarkerFile{$chr."#".$pos_counter}=[0,0,0,0];
				push(@tempMarkerVal,$chr."#".$pos_counter);
				
				#update counts
				if($baseREF eq "A")
				{
					my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
					$temp[0]+=1;
					$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
				}
				elsif($baseREF eq "C")
				{
					my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
					$temp[1]+=1;
					$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
				}
				elsif($baseREF eq "G")
				{
					my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
					$temp[2]+=1;
					$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
				}
				elsif($baseREF eq "T")
				{
					my @temp=@{$tempMarkerFile{$chr."#".$pos_counter}};
					$temp[3]+=1;
					$tempMarkerFile{$chr."#".$pos_counter}=[@temp];
				}
			}
			++$pos_counter;
			++$read_length_counter;
		}
		++$i;
	}
	#print("\n");
	
	my@tempkeys = keys %tempMarkerFile;
	foreach my $val (@tempkeys)
	{
		my @base_counter_cons=@{$Cons{$val}};
		my @temp_array_base_counter=@{$tempMarkerFile{$val}};
		for(my $i=0;$i<@base_counter_cons;++$i)
		{
			$base_counter_cons[$i]+=$temp_array_base_counter[$i];
		}
		$Cons{$val}=[@base_counter_cons];
	}
	#print $begin_." ".$end_." \n";
	return ($begin_,$end_);
}


sub check_rc
{
	my($chr,$position,$read_length)= @_;
	my $hit=0;
	my $hash_entry="";
	for(my $i = $position; $i < ($position+$read_length);++$i)
	{
		#we have a hit! now we have to update the snp_values! here we have to use the marker information
		if(exists $restrictions_site_rs_cons{$chr." ".$i} )
		{
			$hash_entry = $restrictions_site_rs_cons{$chr." ".$i};
			$hit=1;
			last;
		}
	}
	return($hit,$hash_entry);
}



sub readMarkers
{
	my $input_file = new IO::File($marker, "r") or die "could not open $marker: $!\n";
	my $counter = 0;
	while (my $line = $input_file->getline) 
	{
        my @a = split " ", $line;
		my $chr = $a[1];
		$chr=~ s/\D//g; #just add this, because in sorghum the chr are named like  this chr_*
		my $id = $chr."#".$a[2];
		$POS{$id} = [$a[3],$a[4],$a[5],$a[6],$a[7]];
		$Cons{$id}=[0,0,0,0];
		$POS_ref{$id} = $counter;
		++$counter;
	}
}

sub main
{

	my $counter=0;
	my $internal_counter=1;
	my $input_file = new IO::File($mapDotList, "r") or die "could not open $mapDotList: $!\n";
	my $z="";
	if($GZIP==0)
	{
		$z=$input_file;
	}
	else
	{
		$z = new IO::Uncompress::Gunzip $input_file or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
	}
	my $file_writer=IO::File->new();
	$file_writer->open(">".$output_file_name) or die "Could not open ".$output_file_name." $!\n" ;
	
	my $file_writer2=IO::File->new();
	$file_writer2->open(">".$pattern.".txt") or die "Could not open ".$pattern.".txt: $!\n" ;
	my $counter_internal=0;
	
	my %read_pairs=();
	
	while (my $line = $z->getline)
	{		
		chomp($line);
		my @a=split " ",$line;
		if($a[6]!=1  || $a[7]<50)
		{
			next;
		}
		else
		{
			#now we have to check what kind of read we have here, p2 or p1 and if the read hit a rc site!
			
			my $chr = $a[0];
			my $position= $a[1];
			my $pe_flag = $a[9];
			my $read_length = $a[7];
			
			#consider only paired end reads
			if($pe_flag == 3 || $pe_flag == 6)
			{
				if(exists $read_pairs{$a[3]})
				{
					#check distance
					my @temp_array_ = split " ",$read_pairs{$a[3]};
					my $distance=-1;
					if($a[0] == $temp_array_[0])
					{
						$distance = $a[1]-($temp_array_[1]+$temp_array_[7]);
						if($distance <= 420 )
						{
							#distance is ok now we have to check if we have at least on read which mapped towards the rs site!
							my $hit_value=0;
							my $hash_value="";
							if($pe_flag == 3)
							{
								($hit_value,$hash_value) = check_rc($a[0],$position,$read_length);							
							}
							else
							{
								($hit_value,$hash_value) = check_rc($temp_array_[0],$temp_array_[1],$temp_array_[7]);
							}
							#yep we have a valid read pair, now we have to update the snp counter
							if($hit_value==1)
							{
								#update now the nsp
								if($pe_flag == 3)
								{
									update_snp_counter(\@a,\@temp_array_,$hash_value);
								}
								else
								{
									update_snp_counter(\@temp_array_,\@a,$hash_value);
								}
							}								
						}
					}
					delete $read_pairs{$a[3]};
				}
				else
				{
					$read_pairs{$a[3]} = $line;
				}
			}
			else
			{
				next;
			}
		}
	}
	#now sort and print consensus file;
	my @sorted_keys= sort (SortByBeginPos keys %Cons);
	@sorted_keys=sort SortByChr @sorted_keys;
	
	## ALTERNATIVE print allele_count.txt file
	$file_writer->print("#plant_id:$name\n");
    foreach my $val (@sorted_keys)
    {
    	my@a=split "#",$val;
        $file_writer->print($a[0]."\t".$a[1]."\t".(${$POS{$val}}[0]."\t"));
            
		my $ap1 = ${$POS{$val}}[1];
		my $ap2 = ${$POS{$val}}[2];
		my $ap3 = ${$POS{$val}}[3];
		my $ap4 = ${$POS{$val}}[4];

		my $read_number=0;
		for(my $i=0; $i<4;++$i)
		{
			#$file_writer->print(${$Cons{$val}}[$i]."\t"); #A C G T
			$read_number+=${$Cons{$val}}[$i];
		}
		
		$file_writer->print($read_number."\t");
		
		for(my $i=0; $i<4;++$i)
		{
			$file_writer->print(${$Cons{$val}}[$i]."\t"); #A C G T
		}

		for(my $i=1; $i<5;++$i)
		{
			if($i!=4)
			{
				$file_writer->print(${$POS{$val}}[$i]." "); #p1 p2 p3 p4	
			}
			else
			{
				$file_writer->print(${$POS{$val}}[$i]."\n"); #p1 p2 p3 p4
			}	
		}
	}
	
	
	#now we have to print the new file for the basecaller!
	
	
	my @sorted_keys_marker = sort (SortByBeginPos2 keys %restrictions_site_rs_cons_map);
	@sorted_keys_marker = sort SortByChr2 @sorted_keys_marker;
	
	$file_writer->open(">".$output_file_name."_new_format.txt") or die "Could not open "."quality_stat.txt: $!\n";
	$file_writer->print("#plant_id:$name\n");
	foreach my $val (@sorted_keys_marker)
	{
		#get interval
		my $interval_b = $restrictions_site_rs_cons_map{$val}[0];
		if($interval_b eq "A")
		{
			$file_writer->print("% ".$val." N "."\n");
			$file_writer->print("\n");
			next;
		}
		my $interval_e = $restrictions_site_rs_cons_map{$val}[1];
		
		my @a_internal_arry = split " ",$val;
		
		#now get index i and end
		#print " ".$a_internal_arry[0]."#".$interval_b." ... ".$POS_ref{$a_internal_arry[0]."#".$interval_e}."\n";
		#print "val ".$val." begin ".$POS_ref{$a_internal_arry[0]."#".$interval_b}." end ".$POS_ref{$a_internal_arry[0]."#".$interval_e}."\n";
		$file_writer->print("% ".$val." C "."\n");
		for(my $i= $POS_ref{$a_internal_arry[0]."#".$interval_b}; $i<= $POS_ref{$a_internal_arry[0]."#".$interval_e};++$i)
		{
			my $val_internal= $sorted_keys[$i];
			my @a = split "#",$val_internal;
	        $file_writer->print($a[0]."\t".$a[1]."\t".(${$POS{$val_internal}}[0]."\t"));
	            
			my $ap1 = ${$POS{$val_internal}}[1];
			my $ap2 = ${$POS{$val_internal}}[2];
			my $ap3 = ${$POS{$val_internal}}[3];
			my $ap4 = ${$POS{$val_internal}}[4];
	
			my $read_number=0;
			for(my $i=0; $i<4;++$i)
			{
				#$file_writer->print(${$Cons{$val}}[$i]."\t"); #A C G T
				$read_number+=${$Cons{$val_internal}}[$i];
			}
			
			$file_writer->print($read_number."\t");
			
			for(my $i=0; $i<4;++$i)
			{
				$file_writer->print(${$Cons{$val_internal}}[$i]."\t"); #A C G T
			}
	
			for(my $i=1; $i<5;++$i)
			{
				if($i!=4)
				{
					$file_writer->print(${$POS{$val_internal}}[$i]." "); #p1 p2 p3 p4	
				}
				else
				{
					$file_writer->print(${$POS{$val_internal}}[$i]."\n"); #p1 p2 p3 p4
				}	
			}
		}
		$file_writer->print("\n");
	}
		
	#write 
	$file_writer->open(">"."quality_stat.txt") or die "Could not open "."quality_stat.txt: $!\n";
	
	@sorted_keys = keys %quality_base_value;
	
	foreach my $val (@sorted_keys)
	{
		#convert_ascii to int	
		$file_writer->print((ord($val)-33)." ".$quality_base_value{$val}."\n");	
	}
}

sub maxima
{
	my ($a,$b)=@_;
	if($a>$b){return $a;}
	{
		return $b;
	}
}


sub SortByBeginPos2
{
	my @a_array = split " ",$a;
	my @b_array = split " ",$b;
	return $a_array[1]<=>$b_array[1];
}

sub SortByChr2
{
	my @a_array = split " ",$a;
	my @b_array = split " ",$b;
	
	return $a_array[0]<=>$b_array[0];
}


sub SortByBeginPos
{
	my @a_array = split "#",$a;
	my @b_array = split "#",$b;
	return $a_array[1]<=>$b_array[1];
}

sub SortByChr
{
	my @a_array = split "#",$a;
	my @b_array = split "#",$b;
	
	return $a_array[0]<=>$b_array[0];
}

sub addQualityValue
{
	my ($q_string, $counter)=@_;
	if(exists $quality_base_value{substr($q_string,$counter,1)})
	{
		$quality_base_value{substr($q_string,$counter,1)}+=1;
	}
	else
	{
		$quality_base_value{substr($q_string,$counter,1)}=1;
	}	
}
