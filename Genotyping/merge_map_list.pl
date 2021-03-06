use strict;
use warnings;
use Class::Struct;
use Getopt::Std;
use IO::File;
use Time::HiRes;
use Cwd;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
my % options=();
use IO::Handle;
autoflush STDOUT,1;

my$opt_string='hb:';

getopts("$opt_string",\%options) or usage();
usage()if $options{h};

my $stat_barcode_file="";

init();

main();



sub usage
{
	print "perl merge_map.pl -b barcode\n";
}

sub init
{
	if($options{b})
	{
		$stat_barcode_file=$options{b};
		
	}
	else
	{
		die usage();
	}
}



sub main
{
	my $input_file = new IO::File($stat_barcode_file, "r") or die "could not open $stat_barcode_file: $!\n";
	my $cwd = `pwd`;
	chomp($cwd);
		
	while(my $l = $input_file->getline)
	{
		chomp($l);
		my @a=split " ",$l;
	
		if($a[1] >10000)
		{
			my $to_call_string="bsub shore merge -m $a[0]/ -o ../../AlignmentFolder_$a[0] -p ";
			system($to_call_string);
		}
					
	}
	
}
