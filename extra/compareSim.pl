#!/usr/bin/env perl
#
#
#
use strict;
use warnings;
#use 5.010;


my @words=();
my @ref=();
my @prob=();
my @idx=();
my $i=0;
my $randNr=0;
my $prevProb=0;
my $header="";

srand(42);

while (<ARGV>) {
	chop;
	if (not /#/) {
		@words = split(/ /);
		$header = $words[0];
		@words = @words[ 1 .. $#words ];
		@ref = ();
		@prob = ();

		foreach (@words) {
			#print $_."\n";
			if (not /e/) {
				#print $_."\n";
				push @ref,$_;
			} else {
				push @prob,$_;
			}
		}

		@idx = reverse sort { $prob[$a] <=> $prob[$b] } 0 .. $#prob;

		#print @idx;
		
		#print $prob[0];
		#print "\n";

		@ref = @ref[@idx];
		@prob = @prob[@idx];

		#print $_." " foreach @prob;
		#print "\n";
		#print @ref;
		####@prob = @ref[@idx];

		#if ($prob[0]>0.9) {
		
		$i=0;
		$randNr=rand;
		$prevProb=0;
		#print $randNr."\n";
		#while (1==1) {
			#if ($randNr<=$prevProb+$prob[$i]) {
				#last;
				#}
			#	$prevProb=$prevProb+$prob[$i];
			#	$i++;
		#}
		print $header." ".$ref[$i]."\n";
			#}
	}
}
