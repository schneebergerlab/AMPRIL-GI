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

my $output_file="output.txt";
my $file_snp="";
my $dir_input_location="";
my $file_break_name="";

my $BEGIN=5;
my $END=210;

my $GBSFolder ="";
my %hash_snp = ();



my$opt_string='ho:d:g:m:n:';
getopts("$opt_string",\%options) or usage();
usage()if $options{h};

init();
read_snp_marker();
main();

sub usage
{
	my $usage = " usage:perl create_ini_matrix.pl -d dir -m snp_marker file (chr pos p1 p2 p3 p4) -g gbsFolderName -n BREAKFILENAME -o output\n";
	print $usage;
	exit();
}

sub init
{
	if($options{d} && $options{m} && $options{g} && $options{n})
	{
		$dir_input_location=$options{d};
		$file_snp = $options{m};
		$GBSFolder = $options{g};
		$file_break_name = $options{n};
	}
	else
	{
		die usage();	
	}
	
	if($options{o})
	{
		$output_file=$options{o};
	}
}

sub read_snp_marker
{
	my $input_file = new IO::File($file_snp, "r") or die "could not open $file_snp: $!\n";
	
	while(my $line=$input_file->getline)
	{
		chomp($line);
		my @a = split " ",$line;
		
		if(!exists $hash_snp{$a[0]." ".$a[1]})
		{
			print $a[0]." ".$a[1]."\n";
			my @temp = ($a[3],$a[4],$a[5],$a[6]);
			$hash_snp{$a[0]." ".$a[1]}= [@temp];
		}		
	}		
}





sub main
{
	opendir(IMD, $dir_input_location) || die("Cannot open directory");
	my @thefiles= readdir(IMD);
	closedir(IMD);
			
	@thefiles = sort sortByAlph @thefiles;
		
	my $cwd = `pwd`;
	chomp($cwd);
	my %hash_marker_indi=();
	
	my @sorted_values = sort sortByPosition keys %hash_snp; 
	@sorted_values = sort sortByChr @sorted_values; 	
	my @name_array =();
	foreach my $val (@thefiles)
	{
		chdir($cwd);
		if($val =~ m/AlignmentFolder_sample_RAD/)
		{
			print "$val \n";
			#visit each indi get the number of found reads
			my $indi = $val."/$GBSFolder/$file_break_name";
			if(!-e $indi)
			{
				next;
			}
			my $input_file = new IO::File($indi, "r") or die "could not open $indi: $!\n";
			my %hash_genotype_indi_array = ();
			my $counter=0;
			while(my $line = $input_file->getline)
			{
				chomp($line);
				my @a_gento = split " ",$line;
				if($counter==0)
				{
					push(@name_array,$a_gento[0]);
					++$counter;
				}
				
				if(!exists $hash_genotype_indi_array{$a_gento[1]})
				{
					my @temp=($a_gento[2],$a_gento[3]);
					$hash_genotype_indi_array{$a_gento[1]}=[@temp];
				}
				else
				{
					push(@{$hash_genotype_indi_array{$a_gento[1]}},$a_gento[2],$a_gento[3]);
				}				
			}
			#now we have to compare and save the result
			
			my $currentChr=0;
			my @current_genotype_array=0;
			my $genotype_array_counter=0;
			
			foreach my $val (@sorted_values)
			{
				my @val_array = split " ",$val;
				
				if($currentChr eq "" || $currentChr < $val_array[0])
				{
					#get the correct array
					@current_genotype_array = @{$hash_genotype_indi_array{$val_array[0]}};
					#set counter to zero
					$counter=0;
					$currentChr = $val_array[0];
				}
				#find the next big entry or chr change
				if($val_array[1] < $current_genotype_array[$counter])
				{
					my $snp="";
					if($counter == 0)
					{
						 $snp = getAllele($current_genotype_array[$counter+1],\@{$hash_snp{$val}});

					}
					else
					{
						 $snp = getAllele($current_genotype_array[$counter-1],\@{$hash_snp{$val}});
					}
					
					if(!exists $hash_marker_indi{$val})
					{
						my @temp=($snp);
						$hash_marker_indi{$val}=[@temp];														
					}
					else
					{
						push(@{$hash_marker_indi{$val}},$snp);
					}	
				}
				else
				{
					while($counter < @current_genotype_array)
					{	
						if($val_array[1] < $current_genotype_array[$counter])
						{
							last;
						}
						++$counter;
						++$counter;
					}
					my $snp="";
					$snp = getAllele($current_genotype_array[$counter-1],\@{$hash_snp{$val}});
					if(!exists $hash_marker_indi{$val})
					{
						my @temp=($snp);
						$hash_marker_indi{$val}=[@temp];														
					}
					else
					{
						push(@{$hash_marker_indi{$val}},$snp);
					}
				}			
			}		
		}
	}
	
	chdir($cwd);
	my $file_writer=IO::File->new();
	$file_writer->open(">".$output_file) or die "Could not open ".$output_file." $!\n" ;
	
	#first print the name_matrix
	$file_writer->print(" "."\t"." ");
	for(my $i=0; $i < @name_array;++$i)
	{
		$file_writer->print($name_array[$i]."\t");
	}
	$file_writer->print("\n");
	
	foreach my $val (@sorted_values)
	{
		my @snp_values = @{$hash_marker_indi{$val}};
		for(my $i=0; $i < @snp_values;++$i)
		{
			$file_writer->print($snp_values[$i]."\t");
		}
			$file_writer->print("\n");		
	}
			
}

sub getAllele
{
	
	# 0 a 1 b 2 c 3 l
	my ($basecall,$snp)=@_;
	
	#print @$snp ."\n";
	
	if($basecall eq "AA")
	{
		return (@$snp[0]."".@$snp[0]);
	}
	elsif($basecall eq "BB")
	{
		return (@$snp[1]."".@$snp[1]);
	}
	elsif($basecall eq "CC")
	{
		return( @$snp[2]."".@$snp[2]);
	}
	elsif($basecall eq "LL")
	{
		return (@$snp[3]."".@$snp[3]);
	}
	elsif($basecall eq "AC")
	{
		return (@$snp[0]."".@$snp[2]);
	}
	elsif($basecall eq "AL")
	{
		return (@$snp[0]."".@$snp[3]);
	}
	elsif($basecall eq "BC")
	{
		return (@$snp[1]."".@$snp[2]);
	}
	elsif($basecall eq "BL")
	{
		return (@$snp[1]."".@$snp[3]);
	}	
}


sub sortByAlph
{
	return($a cmp $b);
}

sub sortByPosition
{
	my @a_array = split " ",$a;
	my @b_array = split " ",$b;
	return($a_array[1] cmp $b_array[1]);
}

sub sortByChr
{
	my @a_array = split " ",$a;
	my @b_array = split " ",$b;
	return($a_array[0] cmp $b_array[0]);
}
