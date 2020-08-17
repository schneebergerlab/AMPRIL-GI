#! /usr/bin/perl
use strict;
use warnings;
use Class::Struct;
use Getopt::Std;
use IO::File;
use Time::HiRes;
use Cwd;
my % options=();
use IO::Handle;
autoflush STDOUT,1;
my $usage = "$0\n";

###################################
## This needs to be adjusted:

# This builds up the folders in which the consensus summary files are
my $START = 5;
my $END = 220;


# Is the consensus summary file zipped?
my $GZIP = 0;

my $ACOUNTFILE = "allele_count.population_corrected.txt";

my $outfolder = "GBSv1";

my $BSUB_GROUP = "1";

###################################
## This does the job:


my %cross_scheme=();
$cross_scheme{"AxB"}="AxB";
$cross_scheme{"BxA"}="AxB";
$cross_scheme{"ExF"}="ExF";
$cross_scheme{"FxE"}="ExF";
$cross_scheme{"AxC"}="AxC";
$cross_scheme{"CxA"}="AxC";
$cross_scheme{"ExG"}="ExG";
$cross_scheme{"GxE"}="ExG";
$cross_scheme{"AxD"}="AxD";
$cross_scheme{"DxA"}="AxD";
$cross_scheme{"ExH"}="ExH";
$cross_scheme{"HxE"}="ExH";
$cross_scheme{"BxC"}="BxC";
$cross_scheme{"CxB"}="BxC";
$cross_scheme{"FxG"}="FxG";
$cross_scheme{"GxF"}="FxG";
$cross_scheme{"BxD"}="BxD";
$cross_scheme{"DxB"}="BxD";
$cross_scheme{"BxH"}="BxH";
$cross_scheme{"HxB"}="BxH";
$cross_scheme{"CxD"}="CxD";
$cross_scheme{"DxC"}="CxD";
$cross_scheme{"GxH"}="GxH";
$cross_scheme{"HxG"}="GxH";
$cross_scheme{"FxH"}="FxH";
$cross_scheme{"HxF"}="FxH";

#read connection
my $barcode_connection_file="barcode_key.txt";
my %barcode_connetion=();
my $input_file = new IO::File($barcode_connection_file, "r") or die "could not open $barcode_connection_file: $!\n";

while(my $line = $input_file->getline)
{
	chomp($line);
	my @a= split " ",$line;
#	print $a[2]."\n";
	$barcode_connetion{$a[2]} = $cross_scheme{$a[0]};
}
#die;
my $cwd = `pwd`;
chomp($cwd);
print STDERR "Current working dir: $cwd\n";

my $job = 0;

for (my $i = $START; $i <= $END; $i++) {
#foreach my $i (@SAMPLES) {

	# build up folder name:
	my $folder = "AlignmentFolder_sample_RAD-F".$i;
	if(! -e $folder)
	{
		next;
	}
	print STDERR $i, "\t", $folder, "\n";
	chdir($folder);

	if (!-e $outfolder) {
		mkdir($outfolder);
	}
	my $MARKERFILE="/projects/dep_coupland/grp_nordstrom/sequencing/Athal/genome/ampril/Parents/MarkerLists/".$barcode_connetion{$i}."/marker_list_500bp_window";
	my $RESTRICTIONSCUTTINGSITEFILE="/projects/dep_coupland/grp_nordstrom/sequencing/Athal/genome/ampril/Parents/MarkerLists/Restrictionsites/restrictionsites_filtered";
	if (-e "map.list.gz") 
	{
		if ($job == 0) {
			system("bsub -J $BSUB_GROUP.$job \"gunzip map.list.gz\"");
		}
		else
		{
			system("bsub -J $BSUB_GROUP.$job -w 'ended($BSUB_GROUP.".($job-1).")' \"gunzip map.list.gz\"");
			$job++;
		}
	}
	elsif(-e "map.list")
	{
		system("bsub -J $BSUB_GROUP.$job \"perl ~/genotyping/gbs_ampril/allele_count_from_map_list2.pl -n $i -m map.list -p $MARKERFILE -r $RESTRICTIONSCUTTINGSITEFILE -o $outfolder/$ACOUNTFILE "."\"");
		++$job;		
	}
	else {
		die("Could not find map.list for sample_$i\n");
	}
	
	chdir($cwd);
}

