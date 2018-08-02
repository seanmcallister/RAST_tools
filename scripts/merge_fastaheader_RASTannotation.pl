#!/usr/bin/perl -w
#use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
##Take amino acid fasta files from RAST (format: fig|6666666.271858.peg.1).
##Take RAST tab delimited feature file (.txt).
##Convert the header to include: original protein id, RAST annotation, subsystem info, and heme count.
##Use "JJJ" as a field delimiter for later (because lots of scripts will convert all illegal characters to "_").
#
#File organization:
##Dump all feature tables into a single folder. Name of file will match the protein header (6666666.271858.txt).
##All the aa files are in a single folder. Name of file will be sample name (not RAST number).



# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %options=();
getopts("f:a:s:x:o:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-f = path to folder with amino acid fasta files\n";
        print "-x = fasta file suffix\n";
	print "-a = path to folder with RAST tab-delimited txt files\n";
	print "-s = subsystems file from RAST (assets/subsys.txt)\n";
	print "-o = output directory\n";
	print "-h = This help message\n\n";
	print "If you use this script, you should cite the SEED: DOI:10.1093/nar/gki866\n\n";
	die;
    }

my %Sequences;
my %Subsystems;
my %FastaFiles;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN1, "<$options{s}") or die "\n\nNADA $options{s} you FOOL!!!\n\n";
my @subsystems = <IN1>; close(IN1);
foreach my $i (@subsystems)
	{	chomp($i);
		my @subsplit = split('\t', $i);
		if (exists $Subsystems{$subsplit[0]})
			{	$Subsystems{$subsplit[0]}{'Level2'} .= "AND".$subsplit[1];
				$Subsystems{$subsplit[0]}{'Level3'} .= "AND".$subsplit[2];
				$Subsystems{$subsplit[0]}{'Level4'} .= "AND".$subsplit[3];
			}
		else
			{	$Subsystems{$subsplit[0]}{'Level1'} = $subsplit[0];
				$Subsystems{$subsplit[0]}{'Level2'} = $subsplit[1];
				$Subsystems{$subsplit[0]}{'Level3'} = $subsplit[2];
				$Subsystems{$subsplit[0]}{'Level4'} = $subsplit[3];	
			}
	}

#Stores amino acides with protein ids as the dictionary keys.
opendir(DIR, "$options{f}") or die "\n\nNada $options{f} you fool!!\n\n";
my $fastaread = 1;
while (my $file = readdir(DIR))
    {   my $binsuffix = qr/$options{x}/;
        if ($file =~ m/$binsuffix$/)
            {	$FastaFiles{$file}{'Name'} = $file;
		$FastaFiles{$file}{'FileNum'} = $fastaread;
		&FASTAread($file, $fastaread);
                $fastaread += 1;
            }
    }
my $printfastaread = $fastaread - 1;
print "\nTotal number of fasta files: $printfastaread\n";


opendir(DIR2, "$options{a}") or die "\n\nNada $options{a} you fool!!\n\n";
my $rastread = 1;
while (my $file = readdir(DIR2))
    {   if ($file =~ m/\.txt$/)
            {	my $path = $options{a};
		open(IN2, "<".$path."/".$file) or die "\n\nNADA $path"."/"."$file you FOOL!!!\n\n";
		my @RAST = <IN2>; close(IN2); shift(@RAST);
		foreach my $line (@RAST)
			{	my @rastannot = split('\t', $line);
				if (exists $Sequences{$rastannot[1]})
					{	$Sequences{$rastannot[1]}{'start'} = $rastannot[4];
						$Sequences{$rastannot[1]}{'stop'} = $rastannot[5];
						$Sequences{$rastannot[1]}{'strand'} = $rastannot[6];
						$Sequences{$rastannot[1]}{'function'} = $rastannot[7]; #THIS PRIMARY ANNOTATION
						$Sequences{$rastannot[1]}{'aliases'} = $rastannot[8];
						$Sequences{$rastannot[1]}{'figfam'} = $rastannot[9];
						$Sequences{$rastannot[1]}{'evidence_codes'} = $rastannot[10];
						$Sequences{$rastannot[1]}{'nucleotide_sequence'} = $rastannot[11];
					}
			}
                $rastread += 1;
            }
    }
my $printrastread = $rastread - 1;
print "Total number of RAST files: $printrastread\n";


#Calculate heme count
foreach my $i (sort keys %Sequences)
	{	my $aminoseq = $Sequences{$i}{'gappy-ntseq'};
		my $aminoseq2 = $aminoseq;
		$Sequences{$i}{'CXXCH'} = $aminoseq =~ s/C..CH/C..CH/g;
		$Sequences{$i}{'CXXXCH'} = $aminoseq2 =~ s/C...CH/C...CH/g;
	}
foreach my $k (sort keys %Sequences)
	{	unless ($Sequences{$k}{'CXXCH'} >= 1)
		{$Sequences{$k}{'CXXCH'} = 0;}
		unless ($Sequences{$k}{'CXXXCH'} >= 1)
		{$Sequences{$k}{'CXXXCH'} = 0;}
	}


#Match function with subsystem
foreach my $i (sort keys %Sequences)
	{	if (exists $Subsystems{$Sequences{$i}{'function'}})
		  {     $Sequences{$i}{'subsys_lvl2'} .= $Subsystems{$Sequences{$i}{'function'}}{'Level2'};
			$Sequences{$i}{'subsys_lvl3'} .= $Subsystems{$Sequences{$i}{'function'}}{'Level3'};
			$Sequences{$i}{'subsys_lvl4'} .= $Subsystems{$Sequences{$i}{'function'}}{'Level4'};
		  }
	}


#Output
system("mkdir $options{o}");
my $outpath = $options{o};
foreach my $n (sort keys %FastaFiles)
	{	open($n, ">".$outpath."/".$n)
	}

foreach my $fasta (sort keys %Sequences)
	{	my $outhandles = $Sequences{$fasta}{'filename'};
		print $outhandles ">$Sequences{$fasta}{'HEAD'}JJJ$Sequences{$fasta}{'function'}JJJ$Sequences{$fasta}{'subsys_lvl2'}JJJ$Sequences{$fasta}{'subsys_lvl3'}JJJ$Sequences{$fasta}{'subsys_lvl4'}JJJhemecount$Sequences{$fasta}{'CXXCH'}_$Sequences{$fasta}{'CXXXCH'}\n$Sequences{$fasta}{'gappy-ntseq'}\n";
	}

foreach my $n (sort keys %FastaFiles)
	{	close($n);
	}
   
print "\n\n***DONE***\n\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub FASTAread
{	#print "   Reading file . . . \n";
	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
        my $filenumber = $_[1];
        my $path = $options{f};
	open(IN, "<".$path."/".$infile) or die "\n\nNADA $path/$infile you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	#my $unid = $filenumber.10000001;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		$Sequences{$data[0]}{'HEAD'} = $data[0];       # store header
		$Sequences{$data[0]}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		$Sequences{$data[0]}{'SIZE'}    = length($seq);   # store length
		$seq =~ s/\.//;
                $seq =~ s/\-//;
                $Sequences{$data[0]}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $Sequences{$data[0]}{'filenumber'} = $filenumber;
                $Sequences{$data[0]}{'filename'} = $infile;
		$Sequences{$data[0]}{'subsys_lvl2'} = '';
		$Sequences{$data[0]}{'subsys_lvl3'} = '';
		$Sequences{$data[0]}{'subsys_lvl4'} = '';
                #$unid += 1;
	}
	$/="\n";
}
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
