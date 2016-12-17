#!/usr/bin/env perl
#
#
use Data::Dumper qw(Dumper); 
use strict;
use warnings;
use 5.010;
use Getopt::Long qw(:config no_ignore_case);
#use List::MoreUtils qw(uniq);
use Term::Cap;
use POSIX;
#use Storable;

#use List::MoreUtils qw(first_index);

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
open(REFFILE,'<',$reffile) || die "Can't open $reffile : $!\n";

#my @hashArr = @{retrieve($reffile)};

#ref(@hashArr);

#exit;

#$|=1;

my @fields = ();
my $seqDiffPos;
#my $refNr=0;
my $mapPos=0;
my $seq;
my $seqHeader;

my $seq2;
#my @slicenums=();
my @refInd=();
my @curHash;

my @chars;
my @charInd;
my @tmpTot;

my $startScore;
my $curScore;

my $tmp=1;
my @hashArr;
my %tmpHash=();

my $refNr;
my $subNr;

while(my $line=<REFFILE>) {
	chop $line;
	if($line =~ /^>([0-9]+) > ([0-9]+):.*/) {
		$refNr=$1;
		$subNr=$2-1;
		#say $refNr;
		#say $subNr;
	} else {
		#@tmpTot=split(/ /,$line);
		
		#@charInd=map {$tmpTot[$_]} grep {$_ & 1} 1..$#tmpTot;
		#@chars=map {$tmpTot[$_]} grep {$_ & 2} 1..$#tmpTot;

		#$tmpTot[0]=@charInd;
		#$tmpTot[1]=@chars;

		%tmpHash=split(/ /,$line); # interleave indices (keys) and positions (values) onto hash

		$hashArr[$refNr]{$subNr}=();
		for my $key (keys %tmpHash) {
			$hashArr[$refNr]{$subNr}{$key}=$tmpHash{$key};
		}
		#print("\n".$bigArr[$refNr][$subNr]{8749});
		#print Dumper \@bigArr;
		#$tmp++;
		#if ($tmp==3) {exit;}
		#print Dumper \$hashArr[0]{0};
		#exit;
	}
}

#print Dumper \$hashArr[0];

#exit;
#say $hashArr[0][0];

while(my $line=<INFILE>) {
	chop $line;

	if ($line =~ /^@/) {

	} elsif ($line =~ /.*AS:i:([0-9]+)/) {
		$startScore=$1;

		@fields=split(/\t/,$line);
		#print $fields[3]." ";
		
		$refNr=$fields[2]-1;
		$mapPos=$fields[3];
		$seq=$fields[9];
		$seqHeader=$fields[0];

		######
		# get seq2 from cigar string
		#####
		#@curHash=$hashArr[$refNr];
		#print(ref($hashArr[$refNr]));

		#for my $href ( $hashArr[$refNr] ) {
		#	ref($href);
		#     print "{ ";
			#foreach my $x ($href) {
			#	ref($x);
			#}    
			    
			    #for my $role ( keys $href[0] ) {
			    #    print "$role=$href->{$role} ";
			    #}
		#     print "}\n";
		#}

		print $seqHeader." ";
		print $refNr." ";#.$startScore."\n";
		#foreach sub reference in generalized reference with id $refNr:
		
		#print(ref($hashArr[$refNr]));

		#for my $key ( $hashArr[$refNr] ) {
		#for my $i ( 0 .. $#hashArr[$refNr] ) {
			#print(ref($hashArr[$refNr][$i]));
			#foreach my $role (keys %$href) {
			#	print($role=$href->{$role});
			#}
		for my $href ($hashArr[$refNr]) {
			my $bigScore=-10000;
			my $bigKey=0;
		 	for my $key (keys %$href) {
		 		my %curHash=%{$href->{$key}};
			 	#print(Dumper(keys(%curHash)));
			 	#for my $key (keys $hashArr[$refNr]) {
		 		$curScore=$startScore;
		 		#%curHash=$hashArr[$refNr]{$key};

		 		# get changed positions within read frame
		 		@refInd = grep { exists($curHash{$_}) } ($mapPos..($mapPos+99));
		 		#@refInd=keys(%curHash);
		 		#@refInd = @refInd[$mapPos..$tmp2];
		 		
		 		#print $tmp2."\n";
		 		#print scalar @refInd."\n";
		 		#print(Dumper(@refInd));
				#print $curHash{$refInd[0]}." ".substr($seq,$refInd[0]-$mapPos,1)."\n";		 		
	
				# get new alignment score
		 		#for (my $i=0; $i < scalar @refInd; $i++) {
		 		foreach my $ind (@refInd) {
		 			#$seqDiffPos=$mapPos-$refInd[$i];
		 			#if ($seq[$seqDiffPos] ne $curHash{$refInd[$i]}) {
		 			adjustScore(substr($seq,$ind-$mapPos+1,1),$curHash{$ind});
		 			#}
		 		}
		 		#print $key." ".$curScore."\n";
		 		if ($curScore>$bigScore) {
		 			$bigScore=$curScore;
		 			$bigKey=$key;
		 		}
				print "$key $curScore ";
			}
			#print $bigKey;
		}
		print "\n";
	}# elsif ($line =~ /^@/) {
	#}
}

sub adjustScore {
	my $seq=shift;
	my $ref=shift;
	my $penalty=1;

	#print $seq." ".$ref."\n";

	if ($seq eq $ref) { $penalty+=2; }
	elsif ($seq eq '-') { $penalty+=-5; }
	elsif ($seq ne $ref) { $penalty+=-2; }

	$curScore=$curScore+$penalty;
}
