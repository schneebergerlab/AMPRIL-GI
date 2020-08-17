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
my $dir_input_location="";

my$opt_string='ho:d:';
getopts("$opt_string",\%options) or usage();
usage()if $options{h};

init();
main();

sub usage
{
	my $usage = " usage:perl read_number_per_indi.pl -d dir -o output\n";
	print $usage;
	exit();
}

sub init
{
	if($options{d})
	{
		$dir_input_location=$options{d};
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


sub main
{
	opendir(IMD, $dir_input_location) || die("Cannot open directory");
	my @thefiles= readdir(IMD);
	closedir(IMD);
	
	
		
	@thefiles = sort sortByAlph @thefiles;
	
	my %read_counts_per_indi=();
	my %barcode_stat=();
	my $cwd = `pwd`;
	chomp($cwd);
	foreach my $val (@thefiles)
	{
		chdir($cwd);
		if($val =~ m/sample/)
		{
			print "$val \n";
			#visit each indi get the number of found reads
			my $indi = $val."/1/reads_0.fl.gz";
			my $read_counter=0;
			my $input_file = new IO::File($indi, "r") or die "could not open $indi: $!\n";
			my $z = new IO::Uncompress::Gunzip $input_file or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
			
			while(my $line=$z->getline)
			{
				++$read_counter;
			}
			$input_file->close;
			$indi = $val."/2/reads_0.fl.gz";
			$input_file = new IO::File($indi, "r") or die "could not open $indi: $!\n";
			$z = new IO::Uncompress::Gunzip $input_file or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
			while(my $line=$z->getline)
			{
				++$read_counter;
			}
			
			$barcode_stat{$val}=$read_counter;			
			
			if(exists $read_counts_per_indi{$read_counter})
			{
				$read_counts_per_indi{$read_counter}=$read_counts_per_indi{$read_counter}+1;
			}
			else
			{
				$read_counts_per_indi{$read_counter}=1;
			}		
		}
	}
	
	my $file_writer=IO::File->new();
	$file_writer->open(">".$output_file) or die "Could not open ".$output_file." $!\n" ;
	
	my @sorted_values = sort {$a<=>$b} keys %read_counts_per_indi;
	
	foreach my $val (@sorted_values)
	{
		$file_writer->print($val."\t".$read_counts_per_indi{$val}."\n");
	}
	
	$file_writer->close;
	
	$file_writer->open(">".$output_file."_barcode") or die "Could not open ".$output_file."_barcode"." $!\n" ;
	
	@sorted_values = sort {$a cmp $b} keys %barcode_stat;
	
	foreach my $val (@sorted_values)
	{
		$file_writer->print($val."\t".$barcode_stat{$val}."\n");
	}	
		
}

sub sortByAlph
{
	return($a cmp $b);
}
