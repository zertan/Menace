#!/usr/bin/env perl
# generalizeReferences returns a generalized reference from a 
# set of pre aligned LCBs using for instance progressiveMauve
# or mugsy.
#
# input: pre aligend LCBs sorted, with scalar identifiers in
# 		 xmfa format
#
# output: multi-fasta (gaps stripped) 
#
# Daniel Heramansson, 27-09-16

use strict;
use warnings;
use 5.010;
use Getopt::Long qw(:config no_ignore_case);
#use List::MoreUtils qw(uniq);
use Term::Cap;
use POSIX;
use Storable;

my %parserHash=('A'=>'A','C'=>'C','G'=>'G','T'=>'T',
				'AC'=>'M','ACG'=>'V','ACGT'=>'N','ACT'=>'H',
				'AG'=>'R','AGT'=>'D','AT'=>'W', 
				'CG'=>'S','CGT'=>'B','CT'=>'Y',
				'GT'=>'K');

my %rParserHash= reverse %parserHash;
$rParserHash{'-'} .= '';

my %opts;
GetOptions(
	'h|help'        => \$opts{'h'},
	'i|input=s'      => \$opts{'i'},
	'o|output=s'      => \$opts{'o'},
	'r|refs=s'      => \$opts{'r'},
) or die "Error in command line arguments\n";

my $infile = $opts{'i'};
my $outfile = $opts{'o'};
my $reffile = $opts{'r'};


open(INFILE,$infile) || die "Can't open $infile : $!\n";
open(OUTFILE,'>',$outfile) || die "Can't open $outfile : $!\n";
open(REFFILE,'>',$reffile) || die "Can't open $reffile : $!\n";

#$|=1;

my @inseqs=("");
my $seqInd=-1;
my $refNr=-1;
my @headers=();
my @seqDiff=();
my @hashArr=();

#my $seq="";

# anonymous hash
# push @AoH, { husband => "fred", wife => "wilma", daughter => "pebbles" };
# store
# use Storable;
# store \%table, 'file';
# retrieve
# my %host_info = %{retrieve('file')};

while(my $line=<INFILE>) {
	chop $line;
	if($line =~ /^>/) {
		$seqInd++;
		$headers[$seqInd]=$line;
		#print $line;
		#print $seqInd;
		#$inseqs[$seqInd]=$line;
	} elsif ($line =~ /^\#/) {
		#print("cm match");
	} elsif ($line =~ /^=/) {
		$refNr++;
		print OUTFILE ">$refNr\n";
		#say(length($inseqs[0]));
		for (my $i=0; $i < length($inseqs[0]); $i++) {
			getGeneralizedFastaChar($i);
			if (($i+1)%80 == 0) { print OUTFILE "\n"; }
		}
		print OUTFILE "\n\n";
		
		for (my $i=0; $i < scalar @inseqs; $i++) {
			#$hashArr[$refNr][$i]=$seqDiff[$i];
			print REFFILE ">".$refNr." ".$headers[$i]."\n";
			print REFFILE $seqDiff[$i]."\n";
		}

		#store \%table, 'file';
 		#$hashref = retrieve('file');

		$seqInd=-1;
		@inseqs = ("");
		@seqDiff=();
		@headers=();
	} else {
		#print($line);
		#sleep(1);
		$inseqs[$seqInd].=$line;
	}
}

close(INFILE);
close(OUTFILE);
close(REFFILE);

#store \@hashArr, $reffile;
### subs

sub getGeneralizedFastaChar {
	my @chars;


	my $ind=shift;
	my $len=scalar @inseqs;
	$chars[$len-1]="";
	#my $prevCharEqual=0;
	my $first=$rParserHash{substr($inseqs[0],$ind,1)};;

	for (my $i=0; $i < $len; $i++) {
		#print(substr($inseqs[$i],$ind,1));
		$chars[$i]=$rParserHash{substr($inseqs[$i],$ind,1)};
		if ( $first ne $chars[$i]) { $first="False"; }
	}
	if ($first ne "False") { print OUTFILE $first; return; }

	# local $SIG{__WARN__} = sub {
	# 	#die;
	# 	for (my $i=0; $i < $len; $i++) {
	# 		print(substr($inseqs[$i],$ind,1));
	# 	}
	# 	die;
	# };

	#print @chars;

	#return "A";
	@chars=split('',join('',@chars));

	#say scalar @chars;
	#@chars=map(split());
	# eval {
	# 	local $SIG{__WARN__} = sub {
 #    		die;
 #  		};
 #        @chars=split('',join('',@chars));
 #        1;
	# } or do {
 #        my $e = $@;
 #    	print("Something went wrong: $e\n");
	# };
	
	#@chars=map()
	#foreach my $char (@chars) {
        
        #$str=~tr/$key/$rParserHash{$key}/;
    #}
    my $str=join('',sort(uniq(@chars)));
    #$str=~s/(.)(?=.*?\1)//g;
    #$str=sort($str);
	
    #say $str;
	#$str=~tr/-/X/d;
	#return "A";
	#print $ind;
	print OUTFILE $parserHash{$str};
	
	# save difference to $seqDiff
	my $tmp;
	for (my $i=0; $i < $len; $i++) {
		$tmp=substr($inseqs[$i],$ind,1);
		if ($parserHash{$str} ne $tmp) {
		$seqDiff[$i].=$ind." ".$tmp." ";
		#$seqDiff[$i]{$ind} = $tmp;
		}
	}

	#say "";
	#return $parserHash{$str};
	#return;
}

sub uniq {
	my %seen;
	return grep { !$seen{$_}++ } @_;
}
