#! /usr/bin/perl

use strict;
use warnings;
use Class::Struct;
use Getopt::Std;
use IO::File;
use Time::HiRes;
use Cwd;
my $usage = "$0 BSUB_GROUP\n";

###################################
## This needs to be adjusted:

# This builds up the folders in which the consensus summary files are
my $START = 5;
my $END = 48;
# Is the consensus summary file zipped?
my $GZIP = 0;

## Allele count file name prefix

#my $ALLELE_COUNT_PRE = "allele_count.complete";
my $ALLELE_COUNT_PRE = "allele_count.population_corrected.txt_new_format.txt";

my $ALLELE_COUNT_COMPLETE = "allele_count.complete.txt";

# Marker file with the positions for the subsetting
#my $MODELFOLDERBASE = $ENV{HOME}."/genotyping/hmm_profiles";

my $CHRSIZES = "/projects/dep_coupland/grp_nordstrom/data/Athal/TAIR10/chrsizes.chr1-5.txt";

my $OUTPUTFOLDER = "GBSv1";

my $BSUB_GROUP = shift or die $usage;

#my $READ_COUNT_FILE = "reads.raw_count.mapped_count.txt";
#my $MIN_READ_COUNT = 0.025;

#my $BLACK_LIST = "/projects/dep_coupland/grp_nordstrom/sequencing/Athal/genome/recombination/run_134/sample_black_list.txt";

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
	#print $a[2]."\n";
	$barcode_connetion{$a[2]} = $cross_scheme{$a[0]};
}

my $cwd = `pwd`;
chomp($cwd);
print STDERR "Current working dir: $cwd\n";

my $job = 0;

for (my $i = $START; $i <= $END; $i++) {
#foreach my $i (@SAMPLES) {

	# build up folder name:
	my $folder = "AlignmentFolder_sample_RAD-F".$i."/$OUTPUTFOLDER/";
	print STDERR $i, "\t", $folder, "\n";
	chdir($folder);
	my $job = 0;
	my $com = "bsub -q multicore40 -J $BSUB_GROUP.$i.$job -R \"rusage[mem=1000]\" java -Xmx1G -jar /projects/dep_coupland/grp_nordstrom/projects/Genotyping/data/sim_data/population/4parents/jar/base_caller.jar -r $ALLELE_COUNT_PRE -o $ALLELE_COUNT_PRE.base_call.txt -n multi -m";
	system($com);
	++$job;
	
	#hmm
	system("bsub -q short -J $BSUB_GROUP.$i.$job -w 'ended($BSUB_GROUP.$i.".($job-1).")' \"java -jar /projects/dep_coupland/grp_nordstrom/projects/Genotyping/data/sim_data/population/4parents/jar/hmm.jar -r $ALLELE_COUNT_PRE.base_call.txt -o $ALLELE_COUNT_PRE.base_call.hmm.txt -t multi -z /projects/dep_coupland/grp_nordstrom/sequencing/Athal/genome/ampril/MappingPopulation/HMMPROFILE/hmm_profile \"\n"); 
	++$job;
	#conversion
	system("bsub -q short -J $BSUB_GROUP.$i.$job -w 'ended($BSUB_GROUP.$i.".($job-1).")' \"perl ~/genotyping/gbs_ampril/hmm_output_to_marker.pl -n $barcode_connetion{$i} -m $ALLELE_COUNT_PRE -p $ALLELE_COUNT_PRE.base_call.hmm.txt -o $ALLELE_COUNT_PRE.base_call.hmm_marker.txt\"\n");
	#create break file
	++$job;
	system("bsub -q short -J $BSUB_GROUP.$i.$job -w 'ended($BSUB_GROUP.$i.".($job-1).")' \"perl ~/genotyping/gbs_ampril/create_break_file.pl -m$ALLELE_COUNT_PRE.base_call.hmm_marker.txt -c $CHRSIZES -o $ALLELE_COUNT_PRE.base_call.hmm_marker_break.txt\"\n");
	#-m hmm_output_to_marker_file -c chrsize_file -o output_file
	chdir($cwd);
	
	
}
