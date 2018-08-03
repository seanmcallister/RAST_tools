# RAST tools

A collection of tools for converting RAST and MG-RAST outputs for use in other programs. When you call the program use ```-h``` to see the help information.

Tools:

1. ```RAST_to_ANVIO_genetable.pl```
2. ```RAST_to_ANVIO_importannotation.pl```
3. ```RAST_to_GFF.pl```
4. ```merge_fastaheader_RASTannotation.pl```
5. ```mgrast_cov_header.pl```

Tools 1 & 2 are for converting RAST tab-delimited gene files into gene and annotation tables for import into Anvio databases.

Tool 3 converts the RAST tab-delimited gene file into a [GFF-formatted](https://useast.ensembl.org/info/website/upload/gff.html) file for use in bedtools or similar. I use it for identifying the number of reads recruiting to each gene region.

Tool 4 merges information from RAST to make a more informative header per gene. Information includes original protein id, RAST annotation, subsystem info, and heme count. If you use this tool, you should cite the SEED: [DOI:10.1093/nar/gki866](https://academic.oup.com/nar/article/33/17/5691/1067791)

Tool 5 uses the output from the script ```jgi_summarize_bam_contig_depths```, which can be found bundled with [MetaBat](https://bitbucket.org/berkeleylab/metabat) to add coverage information to fasta file headers. This is a requirement for MG-RAST imports.
