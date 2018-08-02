#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Given a directory with binned and unbinned contigs and a directory with RAST (binned) and MG-RAST (unbinned) tab-delimited gene files.
#The assumption is that the contig headers will match the scaffold IDs from the RAST files.
#   1. Test for: Duplicated contigs in the bin directory. Print conflict and die.
#   2. Test for: All contigs found in contig directory have corresponding RAST gene annotations. Print conflict warning.
#   3. Test for: Partition RAST and MG-RAST according to the bins found in the contig directory (even if they are different from the RAST files). Print conflict warning.
#   4. Import RAST files, correctly partitioned by bin or to unbinned.
#   5. Print simple (basic GFF) and annotation enriched GFF files (i.e. annotation, bin or MGrast organism for unbinned, connect to pangenome PC-numbers for Zetas, featureid (peg)).

# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("b:r:o:p:ush", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-r = Folder with RAST tab-delimited gene files with .txt extension. Assumption: RAST (Spreadsheet [tab-separated text format]) and MGRAST files (*_SEED_annotation.function & *_SEED_annotation.organism)\n";
	print "-b = Folder with FINAL bins and unbinned contigs with .fa extension.\n";
        print "-o = Output directory.\n";
        print "-p = Output prefix.\n";
        print "-u = Turn on using unbinned MGRAST info with entire MGRAST pipeline.\n";
        print "-s = Turn on using gene calling for unbinned contigs.\n";
        print "-h = This help message\n\n";
	die;
    }

my %RASTDict;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#1. Test for: Duplicated contigs in the bin directory. Print conflict and die.
print "Testing for duplicated contigs...\n";
my @duplicated_contigs = `cat $options{b}/*.fa | grep ">" | sort | uniq -c | grep -v "1 >"`;
if ($#duplicated_contigs >= 0)
    {   print STDERR "\n\nError: identical contigs found in multiple bins\nOffending contigs:\n";
        foreach my $i (@duplicated_contigs)
            {   chomp($i);
                print STDERR "$duplicated_contigs[$i]\n";
            }
        die "\ndead\n";
    }

#2. Test for: All contigs found in contig directory have corresponding RAST gene annotations. Print conflict warning.
print "Testing for lack of annotation...\n";
my $command1 = "cat $options{b}/*.fa | grep '>' | sed -E 's/^>//' > $options{b}/all_contig_headers_plout.txt";
system("$command1");
my $command2_part1 = "for f in $options{r}/*.txt; ";
my $command2_part2 = 'do tail -n +2 ${f} ';
my $command2_part3 = ">> $options{r}/rast_minus_header.all; done";
my $final_command2 = $command2_part1.$command2_part2.$command2_part3;
system("$final_command2");
my $command3 = "tail -n +2 $options{r}/*.function | cut -f1 | sed -E 's/.+\|(.+_NODE_[0-9]+)_.+/\1/' >> $options{r}/rast_unbinned_header.all";
if ($options{u}) {system("$command3");}
my $command4_part1 = "cat $options{r}/rast_minus_header.all ";
my $command4_part2 = "";
if ($options{u}) {my $command4_part2 = "$options{r}/rast_unbinned_header.all ";}
my $command4_part3 = "| grep 'NODE' | cut -f1 | sort | uniq > $options{b}/all_annotated_headers_plout.txt"; #ASSUMPTION HERE ON CONTIG HEADER NAMES -- consider revising; makes the unannotated count inaccurate
my $final_command4 = $command4_part1.$command4_part2.$command4_part3;
system("$final_command4");
system("cd $options{b}; list_compare.pl -1 all_contig_headers_plout.txt -2 all_annotated_headers_plout.txt 1> /dev/null");
my @missingfrom1 = `cat $options{b}/missingfrom1.txt`;
my @missingfrom2 = `cat $options{b}/missingfrom2.txt`;
if ($#missingfrom1 >= 0)
    {   print STDERR "\n\nError: Gene annotations exist without having a corresponding contig in the bin folder\n";
        print STDERR "Problematic contigs with annotation:\n";
        foreach my $i (@missingfrom1)
            {   chomp($i);
                print STDERR "$i\n";
            }
        die "\ndead\n";
    }
if ($#missingfrom2 >= 0)
    {   system("cat $options{b}/*.fa > $options{b}/allfiles.seanie; grep_ids $options{b}/missingfrom2.txt $options{b}/allfiles.seanie 1> $options{b}/contigs_without_annotation.testlater 2> /dev/null; rm $options{b}/allfiles.seanie");
        system("samtools faidx $options{b}/contigs_without_annotation.testlater");
        print STDERR "\n\nError: contigs are missing gene annotations (not necessarily a problem)\n";
        print STDERR "Potentially problematic contigs have been saved to $options{b}/contigs_without_annotation.testlater\n";
        my $totalproblem = `grep -c ">" $options{b}/contigs_without_annotation.testlater`;
        chomp($totalproblem);
        print STDERR "Summary stats for problematic non-annotated contigs: (n = $totalproblem)\n";
        system("grep '>' $options{b}/contigs_without_annotation.testlater | sed -E 's/^>//' > $options{b}/contigs_without_annotation_headers.testlater");
    if ($options{u} || $options{s}) {    
        my $unbinnedproblem = `grep_ids $options{b}/contigs_without_annotation_headers.testlater $options{b}/*unbinned* 2> /dev/null | grep -c ">"`;
        chomp($unbinnedproblem);
        my $percunbinned = 100 * ($unbinnedproblem / $totalproblem);
        my $rounded = &ROUND($percunbinned,2);
        print STDERR "Percentage of problematic contigs from the unbinned file: $rounded%\n";
        }
        system("cat $options{b}/contigs_without_annotation.testlater.fai | cut -f2 | ~/software/sean_scripts/basic_stats.R");
        print STDERR "\n\n";
    }

#Import RAST files
print "Importing binned RAST annotation files...\n";
open(IN, "<$options{r}/rast_minus_header.all") or die "\n\nFile $options{r}/rast_minus_header.all does not exist.\n\n";
my @data = <IN>; close(IN);
my $unid = 1;
foreach my $line (@data)
	{   chomp($line);
            my @geneinfo = split('\t', $line);
            $RASTDict{$unid}{'contig_id'} = $geneinfo[0];
            $RASTDict{$unid}{'feature_id'} = $geneinfo[1];
            $RASTDict{$unid}{'type'} = $geneinfo[2];
            $RASTDict{$unid}{'strand'} = $geneinfo[6];
            if ($geneinfo[6] eq "+")
                {   $RASTDict{$unid}{'start'} = $geneinfo[4];
                    $RASTDict{$unid}{'stop'} = $geneinfo[5];
                }
            if ($geneinfo[6] eq "-")
                {   $RASTDict{$unid}{'start'} = $geneinfo[5];
                    $RASTDict{$unid}{'stop'} = $geneinfo[4];
                }
            $RASTDict{$unid}{'annotation'} = " ";
            $RASTDict{$unid}{'annotation'} = $geneinfo[7];
            $RASTDict{$unid}{'organism'} = " ";
            $unid += 1;
	}
my $endkey_RAST = $unid - 1;

#Annotate the unbinned file with gene predictions
if ($options{u} || $options{s})
 {  my $startkey_unbinned = $unid;
    print "Running prodigal...\n";
    system("prodigal -i $options{b}/*unbinned* -d $options{b}/UB_genes.fna -p meta -q 1> /dev/null");
    my @ungdata = `grep ">" $options{b}/UB_genes.fna`;
    print "Importing unbinned gene calls from prodigal...\n";
    foreach my $header (@ungdata)
        {   chomp($header);
            $header =~ m/^>((.+_NODE_[0-9]+)_[0-9]+) \# ([0-9]+) \# ([0-9]+) \# ((-|)1) \#.+/;
            $RASTDict{$unid}{'contig_id'} = $2;
            $RASTDict{$unid}{'feature_id'} = $1;
            $RASTDict{$unid}{'type'} = "prodigal_gene";
            if ($5 eq "-1")
                {   $RASTDict{$unid}{'strand'} = "-";   
                }
            elsif ($5 eq "1")
                {   $RASTDict{$unid}{'strand'} = "+";
                }
            $RASTDict{$unid}{'start'} = $3;
            $RASTDict{$unid}{'stop'} = $4;
            $RASTDict{$unid}{'annotation'} = " ";
            $RASTDict{$unid}{'organism'} = "Unbinned";
            $unid += 1;
        }
    my $endkey_unbinned = $unid - 1;

#Remove dictionary entries from binned RAST entries that were really found in the unbinned file
print "Removing duplicated entries from contigs annotated by RAST, but ended up in unbinned.fa...\n";
my %ThrowAway;
foreach my $i (sort keys %RASTDict)
    {   $ThrowAway{$RASTDict{$i}{'contig_id'}}{'RAST_keys'} .= $i."~";
    }
my @toremove;
foreach my $i (sort keys %ThrowAway)
    {   my $keyofint = $ThrowAway{$i}{'RAST_keys'};
        my @keysofint = split("~", $keyofint);
        my $hitunbinned = "FALSE";
        my $hitbinned = "FALSE";
        foreach my $j (@keysofint)
            {   if ($j >= $startkey_unbinned && $j <= $endkey_unbinned)
                    {   $hitunbinned = "TRUE";
                    }
                if ($j >= 1 && $j <= $endkey_RAST)
                    {   $hitbinned = "TRUE";
                    }
                
            }
        if ($hitunbinned eq "TRUE")
            {   foreach my $j (@keysofint)
                    {$RASTDict{$j}{'organism'} = "Unbinned";}
            }
        if ($hitunbinned eq "TRUE" && $hitbinned eq "TRUE")
            {   foreach my $j (@keysofint)
                    {   if ($j >= $startkey_unbinned && $j <= $endkey_unbinned)
                            {   push(@toremove, $j);
                            }
                    }
            }
    }

my @uniq_toremove = uniq @toremove;
my $totalentries = $#uniq_toremove + 1;
print "\tTotal duplicated entries = $totalentries\n";
foreach my $j (@uniq_toremove)
    {   delete $RASTDict{$j};
    }
 }
 
#Import unbinned
print "Making contig ID hash...\n";
my %ContigIds;
foreach my $i (sort keys %RASTDict)
    {   $ContigIds{$RASTDict{$i}{'contig_id'}}{'RAST_keys'} .= $i."~";
        $ContigIds{$RASTDict{$i}{'contig_id'}}{'start'} .= $RASTDict{$i}{'start'}."~";
        $ContigIds{$RASTDict{$i}{'contig_id'}}{'stop'} .= $RASTDict{$i}{'stop'}."~";
        $ContigIds{$RASTDict{$i}{'contig_id'}}{'annotation'} .= $RASTDict{$i}{'annotation'}."~";
        $ContigIds{$RASTDict{$i}{'contig_id'}}{'organism'} .= $RASTDict{$i}{'organism'}."~";
    }

if ($options{u})
{
print "Importing MGRAST annotations to unbinned, where available...\n";
my $functionfilename = ""; $functionfilename = `ls $options{r}/*.function`; chomp($functionfilename);
my $organismfilename = ""; $organismfilename = `ls $options{r}/*.organism`; chomp($organismfilename);
open(IN2, "<$functionfilename") or die "\n\nFile $functionfilename for the unbinned function file does not exist. See -h for help.\n\n";
my @data2 = <IN2>; close(IN2); shift(@data2);
open(IN3, "<$organismfilename") or die "\n\nFile $organismfilename for the unbinned organism file does not exist. See -h for help.\n\n";
my @data3 = <IN3>; close(IN3); shift(@data3);
foreach my $line (0..$#data2)
    {   my @funcinfo = split('\t', $data2[$line]);
        my @orginfo = split('\t', $data3[$line]);
        if ($funcinfo[0] ne $orginfo[0])
            {   print STDERR "The order of the .function file does not match the order of the .organism file.\n\n";
                die "\n\ndead\n\n";
            }
        my $starthead = $funcinfo[0];
        $starthead =~ m/.+\|(.+_NODE_[0-9]+)_\[cov=.+\]_([0-9]+)_([0-9]+)_(.+)\|SEED/;
        my $simplehead = $1;
        my $start = $2;
        my $end = $3;
        my $strand = $4;
    if (defined $ContigIds{$simplehead})
    {
        my $potentialkey = $ContigIds{$simplehead}{'RAST_keys'};
        my @potentialkeys = split("~", $potentialkey);
        my $potentialstart = $ContigIds{$simplehead}{'start'};
        my @potentialstarts = split("~", $potentialstart);
        my $potentialend = $ContigIds{$simplehead}{'stop'};
        my @potentialends = split("~", $potentialend);
        my $potential_annotation = $ContigIds{$simplehead}{'annotation'};
        my @potential_annotations = split("~", $potential_annotation);
        my $potential_organism = $ContigIds{$simplehead}{'organism'};
        my @potential_organisms = split("~", $potential_organism);
        if ($#potentialkeys != $#potentialstarts || $#potentialkeys != $#potentialends ||  $#potentialkeys != $#potential_annotations || $#potentialkeys != $#potential_organisms)
            {   print STDERR "\nSomething is up with contig dictionary.\n\n";}
        foreach my $j (0..$#potentialkeys)
            {   if ($start >= $potentialstarts[$j] - 6 && $start <= $potentialstarts[$j] + 6)
                    {   if ($end >= $potentialends[$j] - 6 && $end <= $potentialends[$j] + 6)
                            {   if ($potential_annotations[$j] eq " " && $potential_organisms[$j] eq "Unbinned")
                                    {   if (defined $funcinfo[12])
                                            {   my $ann = $funcinfo[12];
                                                chomp($ann);
                                                $RASTDict{$potentialkeys[$j]}{'annotation'} = $ann;
                                            }
                                        if (defined $orginfo[12])
                                            {   my $or = $orginfo[12];
                                                chomp($or);
                                                $RASTDict{$potentialkeys[$j]}{'organism'} = "Unbinned_".$or;
                                            }
                                    }
                            }
                    }
            }
    }
    }
}
#3. Test for: Partition RAST and MG-RAST according to the bins found in the contig directory
#(even if they are different from the RAST files). Print conflict warning.
#Replace unbinned "organism" declaraction if found in the bins.
print "Testing for correct partitioning of contigs to the final bin files...\n";
opendir(DIR2, "$options{b}") or die "\n\nFile path $options{b} to the bin folder was not given or does not exist. See -h for help.\n\n";
while (my $file = readdir(DIR2))
    {   if ($file =~ m/\.fa$/ && $file !~ m/unbinned/)
            {   my $path = $options{b}."/".$file;
                my @headers = `cat $path | grep ">" | sed -E 's/^>//'`;
                foreach my $i (@headers)
                    {   chomp($i);
                        if (defined $ContigIds{$i})
                            {   my $potentialkey = $ContigIds{$i}{'RAST_keys'};
                                my @potentialkeys = split("~", $potentialkey);
                                my $potential_organism = $ContigIds{$i}{'organism'};
                                my @potential_organisms = split("~", $potential_organism);
                                if ($#potentialkeys != $#potential_organisms)
                                    {   print STDERR "\nSomething is up with contig dictionary.\n\n";}
                                foreach my $e (0..$#potentialkeys)
                                    {   my $simpfilename = $file;
                                        $simpfilename =~ s/\.fa$//;
                                        if ($potential_organisms[$e] =~ m/^Unbinned_/)
                                            {   print STDERR "Annotation for contig $i was found in unbinned MGRAST file, though contig is binned in $file.\n";
                                            }
                                        if ($potential_organisms[$e] ne " " && $potential_organisms[$e] !~ m/^Unbinned/)
                                            {   print STDERR "Attempting to assign organism $simpfilename to contig $i, however it has already been assigned to bin $potential_organisms[$e]\n";
                                            }
                                        if ($potential_organisms[$e] eq " ")
                                            {   $RASTDict{$potentialkeys[$e]}{'organism'} = $simpfilename;
                                            }
                                    }
                            }
                    }
            }
    }

#opendir(DIR3, "$options{b}") or die "\n\nFile path $options{b} to the bin folder was not given or does not exist. See -h for help.\n\n";
#while (my $file = readdir(DIR3))
#    {   if ($file =~ m/\.fa$/ && $file =~ m/unbinned/)
#            {   my $path = $options{b}."/".$file;
#                my @headers = `cat $path | grep ">" | sed -E 's/^>//'`;
#                foreach my $j (sort keys %RASTDict)
#                    {   foreach my $i (@headers)
#                            {   chomp($i);
#                                if ($i eq $RASTDict{$j}{'contig_id'})
#                                    {   if ($RASTDict{$j}{'organism'} eq " ")
#                                        {   $RASTDict{$j}{'organism'} = "Unbinned";
#                                        }
#                                        last;
#                                    }
#                            }
#                    }
#            }
#    }

#Test for lack of organism assignment
print "Testing for lack of bin/organism assignment...\n";
foreach my $i (sort keys %RASTDict)
    {   if ($RASTDict{$i}{'organism'} eq " ")
            {   print STDERR "Contig $RASTDict{$i}{'contig_id'} is not found in any bin file. Strange...\n\n";
                die "\ndead\n";
            }
        if (exists $RASTDict{$i}{'organism'})
            {   if(defined $RASTDict{$i}{'organism'}){}
                else {print STDERR "Contig $RASTDict{$i}{'contig_id'} organism call isn't defined Strange...\n";}    
            }
        else
            {   print STDERR "Contig $RASTDict{$i}{'contig_id'} organism call doesn't exist. Strange...\n";
            }
        if (exists $RASTDict{$i}{'annotation'})
            {   if(defined $RASTDict{$i}{'annotation'}){}
                else {print STDERR "Contig $RASTDict{$i}{'contig_id'} annotation call isn't defined (not even as empty initialization).\n";}
            }
        else
            {   print STDERR "Contig $RASTDict{$i}{'contig_id'} annotation call doesn't exist (not even as empty initialization).\n"; 
            }
    }

#5. Print simple (basic GFF) and annotation enriched GFF files (i.e. annotation, bin or MGrast organism for unbinned).
#Initial thought to add: connect to pangenome Anvio PC-numbers for Zetas, but just do it outside this program!
system("mkdir $options{o}");
system("mkdir $options{o}/w_strandedness");
print "Saving results to GFF format for use with bedtools...\n";
open(OUT1, ">$options{o}"."/"."$options{p}"."_simple.gff");
open(OUT2, ">$options{o}"."/"."$options{p}"."_annotated.gff");
open(OUT3, ">$options{o}"."/w_strandedness/"."$options{p}"."_stranded_simple.gff");
open(OUT4, ">$options{o}"."/w_strandedness/"."$options{p}"."_stranded_annotated.gff");

foreach my $i (sort keys %RASTDict)
    {   print OUT1 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\t$RASTDict{$i}{'feature_id'}\t0\t$RASTDict{$i}{'strand'}\n";
        print OUT2 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\t$RASTDict{$i}{'feature_id'}\t0\t$RASTDict{$i}{'strand'}\t";
        print OUT2 "$RASTDict{$i}{'organism'}\t$RASTDict{$i}{'annotation'}\n";
        if ($RASTDict{$i}{'strand'} eq "+")
            {   print OUT3 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tanti_$RASTDict{$i}{'feature_id'}\t0\t\-\n";
                print OUT3 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tsense_$RASTDict{$i}{'feature_id'}\t0\t$RASTDict{$i}{'strand'}\n";
                print OUT4 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tanti_$RASTDict{$i}{'feature_id'}\t0\t\-\t";
                print OUT4 "$RASTDict{$i}{'organism'}\t$RASTDict{$i}{'annotation'}\n";
                print OUT4 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tsense_$RASTDict{$i}{'feature_id'}\t0\t$RASTDict{$i}{'strand'}\t";
                print OUT4 "$RASTDict{$i}{'organism'}\t$RASTDict{$i}{'annotation'}\n";
            }
        elsif ($RASTDict{$i}{'strand'} eq "-")
            {   print OUT3 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tanti_$RASTDict{$i}{'feature_id'}\t0\t\+\n";
                print OUT3 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tsense_$RASTDict{$i}{'feature_id'}\t0\t$RASTDict{$i}{'strand'}\n";
                print OUT4 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tanti_$RASTDict{$i}{'feature_id'}\t0\t\+\t";
                print OUT4 "$RASTDict{$i}{'organism'}\t$RASTDict{$i}{'annotation'}\n";
                print OUT4 "$RASTDict{$i}{'contig_id'}\t$RASTDict{$i}{'start'}\t$RASTDict{$i}{'stop'}\tsense_$RASTDict{$i}{'feature_id'}\t0\t$RASTDict{$i}{'strand'}\t";
                print OUT4 "$RASTDict{$i}{'organism'}\t$RASTDict{$i}{'annotation'}\n";
            }
    }
close(OUT1); close(OUT2); close(OUT3); close(OUT4);

print "\n\n* * * * * * * DONE * * * * * * *\n\n";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub ROUND
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -