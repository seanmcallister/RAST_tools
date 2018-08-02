#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Modify the RAST tab-delimited gene files for use in import to ANVIO as a gene caller.
#RAST calls this the "Spreadsheet (tab-separated text format)" under downloads.


# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("i:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = infile\n";
	print "-h = This help message\n\n";
	die;
    }

my %GeneDict;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{i}") or die "\n\nFile $options{i} does not exist or was not given. Try -h for the help file.\n\n";
my @data = <IN>; close(IN); shift(@data);
my $unid = 1000001;
foreach my $line (@data)
	{   	chomp($line);
                my @geneinfo = split('\t', $line);
                my $original_contigid = $geneinfo[0];
                $original_contigid =~ s/\./_/g;
                $original_contigid =~ s/\|/_/g;
                $original_contigid =~ s/-/_/g;
                $original_contigid =~ s/_$//;
                $GeneDict{$unid}{'contig_id'} = uc($original_contigid);
                $GeneDict{$unid}{'feature_id'} = $geneinfo[1];
                $GeneDict{$unid}{'type'} = $geneinfo[2];
                my $start_original = $geneinfo[4];
                my $start_zeroindex = $start_original - 1;
                $GeneDict{$unid}{'start'} = $start_zeroindex;
                my $stop_original = $geneinfo[5];
                my $stop_zeroindex = $stop_original - 1;
                $GeneDict{$unid}{'stop'} = $stop_zeroindex;
                $GeneDict{$unid}{'strand'} = $geneinfo[6];
                $GeneDict{$unid}{'annotation'} = $geneinfo[7];
                $unid += 1;
	}

print "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n";

foreach my $i (sort keys %GeneDict)
    {if ($GeneDict{$i}{'type'} eq "peg") {
        my $geneids = $i - 1000000;
        print "$geneids\t$GeneDict{$i}{'contig_id'}\t";
        if ($GeneDict{$i}{'strand'} eq "+")
            {   my $stupid = $GeneDict{$i}{'stop'} + 1;
                print "$GeneDict{$i}{'start'}\t$stupid\tf\t";
                my $number1 = $GeneDict{$i}{'stop'} - $GeneDict{$i}{'start'} + 1;
                my $number2 = $number1 / 3;
                if (int($number2) != $number2)
                    { print "1\t";}
                else {print "0\t";}
            }
        if ($GeneDict{$i}{'strand'} eq "-")
            {   my $stupider = $GeneDict{$i}{'start'} + 1;
                print "$GeneDict{$i}{'stop'}\t$stupider\tr\t";
                my $number1 = $GeneDict{$i}{'start'} - $GeneDict{$i}{'stop'} + 1;
                my $number2 = $number1 / 3;
                if (int($number2) != $number2)
                    { print "1\t";}
                else {print "0\t";}
            }
        print "rast\tv1.0\n";
      }
    }

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
