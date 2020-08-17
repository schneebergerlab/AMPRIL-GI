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

my$opt_string='hm:o:p:n:';

getopts("$opt_string",\%options) or usage();
usage()if $options{h};

my $file_hmm="";
my $file_marker="";
my $file_output="";
my $name = "";
my %hash_marker=();

init();
read_marker_file();
main();

sub usage
{
	print "perl $0 -n name -m marker_file (allele_count.population_corrected.txt_new_format.txt) -p hmm_file -o output_file\n";
}

sub init
{
	if($options{m} && $options{p} && $options{n})
	{
		$file_marker=$options{m};
		$file_hmm=$options{p};
		$name = $options{n};		
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



sub read_marker_file
{
	my $input_file = new IO::File($file_marker, "r") or die "could not open $file_marker: $!\n";
	
	while (my $line = $input_file->getline) 
	{
		if(substr($line,0,1) eq "%")
		{
			chomp($line);
			my @a = split " ",$line;
			if(!exists $hash_marker{$a[1]})
			{
				my @temp=($a[2]);
				$hash_marker{$a[1]}=[@temp];
			}
			else
			{
				push(@{$hash_marker{$a[1]}},$a[2]);			
			}	
		}		
	}
}

sub main
{

	
	my $input_file = new IO::File($file_hmm, "r") or die "could not open $file_hmm: $!\n";
	
	my $file_writer=IO::File->new();
	$file_writer->open(">".$file_output) or die "Could not open ".$file_output." $!\n" ;
		
	my $counter=0;
		
	while (my $line = $input_file->getline)
	{		
		chomp($line);
		if(substr($line,0,1) eq "#")
		{
			my @a=split " ",$line;
			my @a_internal = split ":",$a[1];
			
			my $chr = $a_internal[1];
			my @array_marker = @{$hash_marker{$chr}};
			
			my @basecalled = split " ",$input_file->getline;
			my @predicted = split " ",$input_file->getline;
			my @masked = split " ",$input_file->getline;
			
			#print @array_marker." ".@basecalled." ".@predicted." ".@masked."\n";
			
			for(my $i=0; $i < @array_marker;++$i)
			{
				$file_writer->print($name." ".$chr."\t".$array_marker[$i]."\t".$basecalled[$i]."\t".$predicted[$i]."\t".$masked[$i]."\n");
			}
		}
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
