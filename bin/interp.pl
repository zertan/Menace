#!/usr/bin/env perl
#
#    Copyright (C) 2009 Genome Research Ltd.
#
#    Author: Heng Li <lh3@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

use strict;
use warnings;


###Builds interpolated pileup from SAM file
##@description counts bases between paired ends and piles up single end reads.
##@output, uses a #header for the RNAME and then the number of reads per base
##@author sm8@sanger.ac.uk, Stephen B. Montgomery

##@caveats
##Requires RNAME to have format as per example
##      chromosome:NCBI36:18:1:76117153:1
##      supercontig::NT_113883:1:137703:1
##      clone::AC138827.3:1:149397:1
##Expects simple CIGAR characters, M, I and D
##Expects SAM file to be sorted.
##Expects 0x0010 to mark second read in PE file (as has been the observed case from MAQ output) (important for line 77)

##Verify and read in SAM file
my $sam_file = $ARGV[0];
if(!defined($sam_file)) { die("No sam file defined on arg 1"); }
unless(-f $sam_file) { die("Sam file does not exist: $sam_file"); }
open(SAM, $sam_file) || die("Cannot open sam file");

##Globals
my $current_location = ""; ##Current RNAME being processed
my $current_size = 0; ##Size of sequence region being processed
my $current_position = 1; ##Current base being processed
my $open = 0; ##Number of open reads (PE reads that have not been closed)
my %close = (); ##Hash of closing positions, when the current_position gets to this position it subtracts the
    ##contained value from those open and deletes the indexed position from the hash

my %lenarr=("NC_002950.2", 2343476,
"NC_004116.1", 2160267,
"NC_004663.1", 6260361,
"NC_004757.1", 2812094,
"NC_006347.1", 5277274,
"NC_006449.1", 1796226,
"NC_007482.1", 635328,
"NC_007519.1", 3730232,
"NC_007712.1", 4171146,
"NC_007716.1", 706569,
"NC_007947.1", 2971517,
"NC_007948.1", 5200264,
"NC_008023.1", 1860355,
"NC_008258.1", 4574284,
"NC_008261.1", 3256683,
"NC_008309.1", 2007700,
"NC_008344.1", 2661057,
"NC_008554.1", 4990251,
"NC_008618.1", 2089645,
"NC_008702.1", 4376040,
"NC_008740.1", 4326849,
"NC_008752.1", 5352772,
"NC_008782.1", 4448856,
"NC_008786.1", 5566749,
"NC_009009.1", 2388435,
"NC_009436.1", 4518712,
"NC_009438.1", 4659220,
"NC_009446.1", 1389350,
"NC_009614.1", 5163189,
"NC_009615.1", 4811379,
"NC_009714.1", 1711273,
"NC_009785.1", 2196662,
"NC_009792.1", 4720462,
"NC_009802.1", 2052007,
"NC_010001.1", 4847594,
"NC_010002.1", 6767514,
"NC_010003.1", 2169548,
"NC_010320.1", 2457259,
"NC_010376.1", 1797577,
"NC_010471.1", 1796284,
"NC_010582.1", 2209198,
"NC_010655.1", 2664102,
"NC_011047.1", 601943,
"NC_011283.1", 5641239,
"NC_011662.2", 4496212,
"NC_011740.1", 4588711,
"NC_011769.1", 4040304,
"NC_011883.1", 2873437,
"NC_011898.1", 4068724,
"NC_011992.1", 3796573,
"NC_011999.1", 2102324,
"NC_012691.1", 3471292,
"NC_012778.1", 2144190,
"NC_012913.1", 2313035,
"NC_013118.1", 1220319,
"NC_013192.1", 2465610,
"NC_013203.1", 1543805,
"NC_013204.1", 3632260,
"NC_013316.1", 4191339,
"NC_013504.1", 1755993,
"NC_013508.1", 3760463,
"NC_013515.1", 1662578,
"NC_013520.1", 2132142,
"NC_013522.1", 1848474,
"NC_013714.1", 2636367,
"NC_013715.1", 2264603,
"NC_013721.1", 1617545,
"NC_013740.1", 2329769,
"NC_013798.1", 2350911,
"NC_013850.1", 5458505,
"NC_013853.1", 2146611,
"NC_013895.2", 1809746,
"NC_014033.1", 3619559,
"NC_014106.1", 2043161,
"NC_014246.1", 2146480,
"NC_014363.1", 2051896,
"NC_014365.1", 3655731,
"NC_014370.1", 1796408,
"NC_014387.1", 3554804,
"NC_014618.1", 4814049,
"NC_014624.2", 4316707,
"NC_014632.1", 2046464,
"NC_014656.1", 2265943,
"NC_014724.1", 2067702,
"NC_014752.1", 2220606,
"NC_014800.1", 704884,
"NC_014828.1", 3008576,
"NC_014833.1", 3685408,
"NC_014844.1", 3629109,
"NC_014933.1", 3998906,
"NC_015138.1", 5482170,
"NC_015160.1", 4392288,
"NC_015164.1", 4242803,
"NC_015214.1", 2078001,
"NC_015291.1", 1958690,
"NC_015311.1", 2937589,
"NC_015385.1", 2731853,
"NC_015389.1", 2115681,
"NC_015422.1", 4995263,
"NC_015437.1", 2568361,
"NC_015501.1", 2186370,
"NC_015563.1", 6685842,
"NC_015600.1", 2100077,
"NC_015677.1", 4070193,
"NC_015678.1", 2153652,
"NC_015737.1", 2835737,
"NC_015738.1", 3123671,
"NC_015740.1", 4547930,
"NC_015856.1", 5186898,
"NC_015873.1", 2474718,
"NC_015875.1", 2190731,
"NC_015964.1", 2086875,
"NC_015968.1", 4812833,
"NC_015975.1", 2066652,
"NC_015977.1", 3592125,
"NC_016048.1", 4410036,
"NC_016052.1", 2562720,
"NC_016077.1", 2487765,
"NC_016445.1", 3031375,
"NC_016513.1", 2136808,
"NC_016610.1", 3405521,
"NC_016616.1", 3806980,
"NC_016628.1", 1621862,
"NC_016630.1", 1931012,
"NC_016749.1", 2130034,
"NC_016818.1", 4861101,
"NC_016822.1", 4988504,
"NC_016826.1", 1988420,
"NC_017068.1", 3003680,
"NC_017280.1", 1616648,
"NC_017451.1", 1932306,
"NC_017481.1", 1828169,
"NC_017511.1", 2154835,
"NC_017517.1", 2287777,
"NC_017568.1", 2572216,
"NC_017595.1", 2210574,
"NC_017617.1", 2038034,
"NC_017625.1", 4630707,
"NC_017845.1", 5164411,
"NC_017856.1", 2697465,
"NC_017861.1", 2119790,
"NC_017866.1", 1938595,
"NC_017910.1", 4158725,
"NC_017999.1", 2223664,
"NC_018011.1", 3734239,
"NC_018106.1", 6097032,
"NC_018221.1", 2987450,
"NC_018607.1", 2765477,
"NC_018691.1", 4928223,
"NC_018708.1", 5196935,
"NC_018712.1", 2151145,
"NC_019968.1", 1450390,
"NC_020125.1", 2166321,
"NC_020127.1", 1457568,
"NC_020181.1", 5419609,
"NC_020294.1", 833125,
"NC_020450.1", 2421471,
"NC_020515.1", 2407846,
"NC_020517.1", 2422684,
"NC_020546.1", 2291643,
"NC_021009.1", 3522704,
"NC_021010.1", 3344951,
"NC_021011.1", 2943413,
"NC_021012.1", 4286292,
"NC_021013.1", 2249085,
"NC_021014.1", 3545606,
"NC_021015.1", 3341681,
"NC_021017.1", 5976145,
"NC_021018.1", 3123007,
"NC_021019.1", 1966750,
"NC_021021.1", 3608022,
"NC_021022.1", 3757491,
"NC_021023.1", 3096657,
"NC_021030.1", 3763317,
"NC_021031.1", 3164379,
"NC_021035.1", 3601020,
"NC_021038.1", 2728333,
"NC_021039.1", 2573208,
"NC_021041.1", 2209938,
"NC_021042.1", 3321367,
"NC_021046.1", 4908759,
"NC_021047.1", 3769775,
"NC_021066.1", 5398151,
"NC_021082.1", 2731870,
"NC_021175.1", 2142100,
"NC_021281.1", 2261267,
"NC_021500.1", 5039027,
"NC_021659.1", 5467306,
"NC_021900.1", 1975547,
"NC_021994.1", 2994661,
"NC_022124.1", 709850,
"NC_022196.1", 2268272,
"NC_022198.1", 2031902,
"NC_022239.1", 2233640,
"NC_022245.1", 1935662,
"NC_022246.1", 1996214,
"NC_022247.1", 5795261,
"NC_022567.1", 2862526,
"NC_022576.1", 2991840,
"NC_022582.1", 2023580,
"NC_022584.1", 1992567,
"NC_022604.1", 2987966,
"NC_022738.1", 5644569,
"NC_023032.1", 4377544,
"NC_023061.1", 6683584,
"NT_187110.1", 11127);

while (my $line = <SAM>) {
    my @tokens = split /\t/, $line;

    if ($current_location ne $tokens[2]) { ##Start a new sequence region
        for (my $i = $current_position; $i <= $current_size; $i++) { ##Close the previous sequence region
            if (defined($close{$i})) {
                $open = $open - $close{$i};
                delete $close{$i};
            }
            print $i . " " . $open . "\n";
        }
        if ($current_location ne "") {
            print "\n";
        }

        ##Initiate a new sequence region
        my @location_tokens = split /\|/, $tokens[2];
		$current_position = 1;
        $current_location = $tokens[2];
        $current_size = $lenarr{$location_tokens[5]}; #$location_tokens[4];
		$open = 0;
        %close = ();
		#print "#" . $tokens[2] . "\n";

        ##Print pileup to just before the first read (will be 0)
        for (my $current_position = 1; $current_position < $tokens[3]; $current_position++) {
            print $current_position . " " . $open . "\n";
			#print $open . "\n";
        }
        $current_position = $tokens[3];

    } else { ##Sequence region already open
        if ($tokens[3] > $current_position) { ##If the new read's position is greater than the current position
                                                ##cycle through to catch up to the current position
            for (my $i = $current_position; $i < $tokens[3]; $i++) {
                if (defined($close{$i})) {
					$open = $open - $close{$i};
                    delete $close{$i};
                }
                print $i . " " . $open . "\n";
            }
            $current_position = $tokens[3];
        }
    }
	#$open++; ##Increment the number of open reads
	
	#print $line . "\n";

	#if (($tokens[1] & 0x0080 || $tokens[1] & 0x0040) && $tokens[1] & 0x0020 && $tokens[1] & 0x0002) { ##if second read of mate pair, add close condition
	#	$open++;
#		my $seq_region_end=$tokens[8]+$tokens[3]-1;
#		if (!defined($close{$seq_region_end + 1})) { $close{$seq_region_end + 1} = 0; }
#		$close{$seq_region_end + 1} = $close{$seq_region_end + 1} + 1;
		#$open--;
        
		#my $parsed_cig = &parseCigar($tokens[5]);
        
		#my $seq_region_end = $tokens[3] + $parsed_cig->{'M'} + $parsed_cig->{'D'} - 1;
		#if (!defined($close{$seq_region_end + 1})) { $close{$seq_region_end + 1} = 0; }
		#$close{$seq_region_end + 1} = $close{$seq_region_end + 1} + 1;
		#print "###".$open-$paired."###";
		#$paired++;
		#print "paired close";
		#foreach my $position (keys %close)
		#{
		#		print "$position\n";
		#}
	if (($tokens[1] & 0x0004) || ($tokens[1] & 0x0008)  || ($tokens[1] & 0x0100)  || ($tokens[1] & 0x0200) || ($tokens[1] & 0x0400) || ($tokens[1] & 0x0800)) {
		#print "cond1:\n";
		#print "$line\n";
		#<STDIN>;
		# nothing
	} elsif (!($tokens[1] & 0x0001) || !($tokens[1] & 0x0002)) { ##if unpaired, add close condition
		#print "cond2:\n";
		#print "$line\n";
		#<STDIN>;
		
		$open++;
		my $parsed_cig = &parseCigar($tokens[5]);
		my $seq_region_end = $tokens[3] + $parsed_cig->{'M'} + $parsed_cig->{'D'} - 1;
		if (!defined($close{$seq_region_end + 1})) { $close{$seq_region_end + 1} = 0; }
		$close{$seq_region_end + 1} = $close{$seq_region_end + 1} + 1;
		#		#print "single\n"
	} elsif (( $tokens[1] & 0x0080) && $tokens[1] & 0x0010 && $tokens[1] & 0x0002 && $tokens[1] & 0x0001) {
		#print "cond3:\n";
		#print "$line\n";
		#<STDIN>;
		#if (abs($tokens[8])<600) {
		#	$open++;
		#	my $seq_region_end=$tokens[8]+$tokens[3]-1;
		#	if (!defined($close{$seq_region_end + 1})) { $close{$seq_region_end + 1} = 0; }
		#	$close{$seq_region_end + 1} = $close{$seq_region_end + 1} + 1;
		#}
	} elsif (( $tokens[1] & 0x0040) && $tokens[1] & 0x0020 && $tokens[1] & 0x0002 && $tokens[1] & 0x0001) {
		if (abs($tokens[8])<600 && $tokens[8]>0 && $tokens[7]>$current_position) {
			#print "cond4:\n";
			#print "$line\n";
			#<STDIN>;
			$open++;
			my $seq_region_end=$tokens[8]+$tokens[3]-1;
			if (!defined($close{$seq_region_end + 1})) { $close{$seq_region_end + 1} = 0; }
			$close{$seq_region_end + 1} = $close{$seq_region_end + 1} + 1;
			#print "\n$seq_region_end\n";
		} #else {
			#print "$line\n";
			#}
	}
		#$paired--;
		#print "paired open"
        #do nothing
		# }
	#my $ll=$open-length(keys(%close));
	#my $size = keys %close;
	#print "xx" . $size . "xx";
	#print "\tdiff: $ll\n";
	#<STDIN>;	
}

#print $current_position
#$current_size = 5277274;
for (my $i = $current_position; $i <= $current_size; $i++) {  ##Finish up the last sequence region
    if (defined($close{$i})) {
        $open = $open - $close{$i};
        delete $close{$i};
    }
    print $i . " " . $open . "\n";
	#print $open . "\n";
}
print "\n";
close(SAM);

#foreach my $position (keys %close)
#{
#	print "$position\n";
#}

exit(0);

##reads and tokenizes simple cigarline
sub parseCigar() {
    my $cigar_line = shift;
    $cigar_line =~ s/([0-9]*[A-Z]{1})/$1\t/g;
    my @cigar_tokens = split /\t/, $cigar_line;
    my %parsed = ('M' => 0,
                  'I' => 0,
                  'D' => 0);
    my @events = ();
    for(my $i = 0; $i < scalar(@cigar_tokens); $i++) {
        if ($cigar_tokens[$i] =~ /([0-9]+)([A-Z]{1})/g) {
            if (!defined($parsed{$2})) { $parsed{$2} = 0; }
            my $nt = $2;
            if ($nt ne "M" && $nt ne "D"  && $nt ne "I") { $nt = "M"; }
            $parsed{$nt} += $1;
            my %event_el = ("t" => $nt,
                            "n" => $1);
            push @events, \%event_el;
        }
    }
    $parsed{'events'} = \@events;
    return \%parsed;
}
