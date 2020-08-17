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

my$opt_string='hm:o:c:';

getopts("$opt_string",\%options) or usage();
usage()if $options{h};

my $file_chr="";
my $file_marker="";
my $file_output="";

my %hash_chr=();
#/projects/dep_coupland/grp_nordstrom/data/Athal/TAIR10/chrsizes.chr1-5.txt
init();
read_chr_size();
main();

sub usage
{
	print "perl $0 -m hmm_output_to_marker_file -c chrsize_file -o output_file\n";
}

sub init
{
	if($options{m} && $options{c})
	{
		$file_marker=$options{m};
		$file_chr=$options{c};		
	}
	else
	{
		die usage();
	}	
	if($options{o})
	{
		$file_output=$options{o};
	}
}

sub read_chr_size
{
	my $input_file = new IO::File($file_chr, "r") or die "could not open $file_chr: $!\n";
	
	while (my $line = $input_file->getline) 
	{
		chomp($line);
		my @a = split " ",$line;
		if(!exists $hash_chr{$a[0]})
		{
			$hash_chr{$a[0]} = $a[1];
		}
	}	
}

sub main
{
	my $input_file = new IO::File($file_marker, "r") or die "could not open $file_marker: $!\n";
	
	my $file_writer=IO::File->new();
	$file_writer->open(">".$file_output) or die "Could not open ".$file_output." $!\n" ;
		
	my $counter=0;
		
	my $prev_genotype ="";
	my $position = "";
	my $chr="";
	my $name ="";
	while (my $line = $input_file->getline)
	{		
		chomp($line);
		my @a=split " ",$line;
		if($name eq "")
		{
			$name = $a[0];
		}
		if($chr eq "" || $chr < $a[1])
		{
			if($chr ne "")
			{
				#print $chr." ".$prev_genotype."\n";
				$file_writer->print($name."\t".$chr."\t".$hash_chr{$chr}."\t".$prev_genotype."\n");
			}
			$chr = $a[1];
			$position = 0;
			$prev_genotype = "";
		}
		
		if($prev_genotype eq "")
		{
			$prev_genotype = $a[4];
		}
		else
		{
			#check if we have the same genotype, if not change it + print it!
			if($prev_genotype ne $a[4])
			{				
				my $middle_distance = int ($position-1+($a[2]-$position)/2);
				$file_writer->print($name."\t".$chr."\t".$middle_distance."\t".$prev_genotype."\n");
				$prev_genotype = $a[4];
				$position = $middle_distance+1;
			}
			else
			{
				$position=$a[2];
			}
		}	
	}
	$file_writer->print($name."\t".$chr."\t".$hash_chr{$chr}."\t".$prev_genotype."\n");
}
