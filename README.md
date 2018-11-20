# Syntenizer3000
Tool for bacterial genome analysis including but not limited to synteny measurement.

Run build.sh in the folder with the source code to build an executable on a UNIX derived OS. GCC with at least C++11 support or equivalent must be installed.

Run without any input parameters to get a list of possible arguments and a short explanation of what they do.

Below is an example of running the program on the example dataset with a full suite of analysis enabled:

./Syntenizer --gffdir=/data/input --outdir=/data/output --genegroups=/data/input/gene_groups.csv --contigcolours=/data/input/ContigColours.csv --syntenizegroups --generatestraincharts --generatesyntenycharts --chromosomeidentifier=chromosome --plasmididentifier=Rh --fragmentidentifier=fragment --genegrouphighlightcolours=/data/input/GeneGroupHighlightColours.csv --ffndir=/data/input --generatefnagroups --classifyfragments --disambiguategroups

In the above example all input files, including *.gff, *.fnn gene groups and others are located in the folder '/data/input', and all files generated by the tool are output into the folder '/data/output'.

If the program has been configured to output circular strain chart data or gene group synteny map data then the next step is to render these using the provided R scripts 'GenerateCharts.R', 'GenerateSyntenyMaps.R' and 'GenerateHeatmap.R'. These scripts need to be manually run by users with a single input parameter, which is the path to the folder where the chart or map data generated by Synternizer has been output, which in the above example would be '/data/output'.

Things to note:

If issues or uncertainties crop up while trying to use the program you can refer to the example dataset provided. A full suit of analysis can be performed using said dataset so it should serve as a guide on how input data should be formatted for the program to work.

Only GFF3 files generated by PROKKA have been tested. Other versions or files generated by other tools may contain variations that might cause incompatibility issues. In that case contact the author for possible solutions or fixes.

Two gene group file formats are supported. The direct gene group file output from ProteinOrtho, plus an internal format detailed below.

Each gene group in the internal format is defined on a space (ASCII 0x20) separated line starting with The group id as the first item. All subsequent items should be a concatenation of the strain id and gene id with the character '|' (ASCII 0x7C) separating them. See gene_groups.csv in the example dataset for an example.

To generate basic circular charts the program only needs properly formatted GFF3 files and accompanying gene group file.
However, the GC3s lane requires the presence of gene sequence data. If these are part of the GFF3 files then they will be parsed automatically. Alternatively they can be provided using gene group *.fna files using the input parameter.

In order for contigs to be coloured according to their type (chromosome/plasmid/fragment) then the contig ID's in the GFF3 files have to contain an string patter that identifies them as such. The program can then be configured to recognize these using the appropriate options as detailed below. For examples of this working correctly use the example input parameters provided above on the example dataset.

When Syntenizer is run with a bare minimum of input parameters, which are the *.GFF file input folder, output folder and gene group file three files will be output:

gene_groups.csv, which contains a Syntenizer formatted version of the the input gene groups. If the gene group input file was a Syntenizer formatted file then this output file should be identical to the input file in terms of conntents. Since gene group highlighting on charts requires referrencing to gene group ID's, and since ProteinOrto formatted gene groups lack these it is necessary to run the ProteinOrtho gene groups through Syntenizer and use the obtained Syntenizer formatted gene group output file as a basis for creating charts if gene groups need to be highlighted.

presence_absence_matrix.csv:
A file that contains a presence/absence matrix based on gene groups. The ';' (ASCII 0x3B) delimited file represents which strains are present "1" or absent "0" in each gene group.

presence_absence_gapit.csv:
A GAPIT formatted presence/absence matrix. This file can be directly loaded and passed to GAPIT for use in GWAS analysis of the given bacterial strains.

gene_relation_matrix.csv:
This file contains a gene relation matrix based on presence/absence of strains in gene groups. This file is needed for rendering a gene relation heatmap using the GenerateHeatmap.R script. Consult the manual for information on how the matrix is calculated.

strain_metrics.csv:
This file contains useful statistics on the input dataset. Below is a short explanation of each column.

Strain ID   -   Name of the given strain.
Contigs     -   The number of contigs in the given strain.
bp          -   The number of determined basepairs in the given strain.
ubp         -   The number of undetermined basepairs in the given strain. 
GC          -   The GC-content ratio [0:1] of the given strain.
CDS         -   The number of coding sequences in the given strain. Note that CDS = Orthologs+Paralogs+Uniques.
Orthologs   -   The number of orthologous genes in the given strain.
Paralogs    -   The number of paralogous genes in the given strain
Uniques     -   The number of unique genes in the given strain. 
tRNA        -   The number of tRNA in the given strain, if present in the GFF file.
rRNA        -   The number of rRNA in the given strain, if present in the GFF file.
tmRNA       -   The number if tmRNA in the given strain, if present in the GFF file.
0-10k       -   The percentage of total strain basepairs located on contigs between 0 to 10kbp in size.
10-250k     -   The percentage of total strain basepairs located on contigs between 10 to 250kbp in size.
250k-500k   -   The percentage of total strain basepairs located on contigs between 250 to 500kbp in size.
500k-750k   -   The percentage of total strain basepairs located on contigs between 500 to 750kbp in size.
750k-1000k  -   The percentage of total strain basepairs located on contigs between 750 to 1mbp in size.
1000-5000k  -   The percentage of total strain basepairs located on contigs between 1 to 5mbp in size.
>5000k      -   The percentage of total strain basepairs located on contigs larger than 5mbp in size.

Below is a list and explanation of all input parameter options offered by Syntenizer.

INPUT OPTIONS

--gffdir=dir
Sets the directory where strain gff files are located. All files with the .gff extension will be loaded. The strain ID is           set to the name of the file (with the the extension removed) i.e. the file calb.gff defines the genome of the calb strain.    If contig size and contig sequence data are present in the gff file then these will be parsed. If they are not present the size of each contig is set to the end of the final coding sequence gene defined on said contig. Alternatively sequence data for each gene can be provided using .ffn files. Consult the demo dataset for examples of file formatting. The function of a gene defined using the 'product=' tag will be parsed if it is present. Consult the demo dataset as for examples.
    
--outdir=dir
Sets the output directory for all files generated by Syntenizer3000. This directory must exist beforehand, and any files already present with similar name and extension as output files will be overwritten.

--genegroups=file
Set the file specifying gene groups. Syntenizer3000 can parse two different ways of encoding this information. The first is the direct output file format from the orthology tool ProteinOrtho. The second is an internal format, i.e. as below:
Group1:	strain1|gene1	strain2|gene1
Each tab delimited line starts with a group ID finished with a colon. Then a list of genes 
in said group each consisting of the ID of the gene strain and ID of the gene in said strain, separated by the | ASCII character (0xC7).

--ffndir=dir
Sets the directory where gene sequence .ffn files are located. All files with the ffn extension are parsed, and the gene sequences found within are attached to their corresponding genes. Sequences for unknown genes are ignored. Sequences that may already have been parsed from .gff will be compared to the ones parsed from .ffn files and result in an error if not an exact match.

--contigcolours=file
Sets the file specifying the contig colours in strain charts. Consult the file colours.csv from the demo dataset for an example of how the files is to be formatted. The first tab delimited line is a header defining the role of each subsequent tab delimited line. Each line contains two items, the first is a string pattern, and the second is a valid R colour i.e red, blue, #C0C0C0 etc. The ID of each contig being rendered is searched for the defined string patterns and if there is a match the colour associated with the pattern will be applied to the contig. If the contig ID matches several of the defined string patterns then the colour of the string pattern nearest the end of the contig colour file will be the one used.

--syntenizegroups
Generate synteny score for each gene group into the file groups_synteny.csv, in addition to other useful information about said gene group. Below is an explanation of each column in the ; delimited file:
Group - ID of the group	
Connectivity - The connectivity of the group as parsed from the ProteinOrtho file.
Genes - The number of genes in group.
Strains - The number of unique strains of the genes in group.
Orthologs - The number of orthologous genes in group.
Paralogs - The number of paralogrous genes in the group
Chromosome Genes - The number of genes residing on chromosomal contigs in the group. This stat only works if a valid chromosome identifier. Consult the demo dataset for an example of this functioning correctly.
Plasmid Genes - The number of genes residing on plasmid contigs in the group. This stat only works if a valid plasmid identifier. Consult the demo dataset for an example of this functioning correctly.
Fragment Genes - The number of genes residing on fragment contigs in the group. This stat only works if a valid fragment identifier. Consult the demo dataset for an example of this functioning correctly.
Protein Products - All the protein products of the genes in the group, separated by the / character. The protein product is parsed from .gff files using the 'product=' tag if present. In an ideal world each gene group should only contain one protein product, and differing protein products could be a sign of a badly constructed group. Consult the demo dataset for an example of this functioning correctly.
Synteny Score SIMPLE - A number between 0 and 40 representing the average SIMPLE synteny score of the group. Consult the manual for details on how this is  calciulated.
Synteny Score SOPHISTICATED - A number between 0 and 40 representing the SOPHISTICATED score of the group. Consult the manual for details on how this is calculated.

--disambiguategroups
Generates file disambiguated_gene_groups.csv with gene groups devoid of paralogs. Consult the manual for details on the disambiguation process.

--generatefnagroups
Generate gene group .fna files which define the sequence of each member of the gene group.

--generatestraincharts
Generates the data files needed for rendering strain charts using the GenerateCharts.R script. Consult the manual for information on what data these charts depict.

--generatesyntenycharts
Generate the data files needed for rendering group synteny charts using the GenerateSyntenyMaps.R script. Consult the manual for information on what data these charts depict.

--chromosomeidentifier=pattern
Sets the character pattern that identifies a chromosomal contig if it is found to be part of its ID. This parameter is needed for contig classification to work. Consult the demo dataset for an example of how this could work.

--plasmididentifier=pattern
Sets the character pattern that identifies a plasmid contig if it is found to be part of its ID. This parameter is needed for contig classification to work. Consult the demo dataset for an example of how this could work.

--fragmentidentifier=pattern
Sets the character pattern that identifies a fragment contig if it is found to be part of its ID. Two number digit characters are assumed to exist right after this pattern in any matched contig ID, which identify the plasmid type, i.e. Rh01. This parameter is needed for contig classification to work. Consult the demo dataset for an example of how this could work.

--syntenizegenepairs=file
Sets the file specifying specific gene pairs to be syntenized into the gene_pairs_synteny.csv output file. Consult the file GenePairs.csv in the demo dataset for an example of formatting.

--genegrouphighlightcolours=file
Sets the file specifying highlight colours for gene group in rendered strain charts. Each ; delimited line defines two items, the first being a gene group ID, and the second a valid R colour said gene group should be highlighted with on rendered strain charts. Consult the file GeneGroupHighlightColours.csv in the demo dataset for an example of formatting.

--classifyfragments
Generate file classified_fragments.csv with fragment classifications based on synteny and gene group information. Below is an explanation of each column in the ; delimited file:
Strain - The ID of the strain this contig belongs to.
Contig - The ID of the contig.
Genes - The number of coding genes on this contig.
Chromosome - A number between 0 and 1 describing the fraction of chromosomal genes in the gene groups of the genes on this contig.
Plasmid - A number between 0 and 1 describing the fraction of plasmid genes in the gene groups of the genes on this contig.
Chromosomal Synteny - A number between 0 and 40 describing the average synteny between chromosomal genes in the gene groups of each gene on this contig.
Plasmid Synteny - A number between 0 and 40 describing the average synteny between plasmid genes in the gene groups of each gene on this contig.
Coverage - A number between -1 and 1 calculated as "chromosome - plasmid". Consult the manual for an explanation of this parameter.
Syneny - A number between -1 and 1 calculated by "(Chromosomal Synteny / 40) - (Plasmid Synteny / 40)". Consult the manual for an explanation of this parameter.
Magnitude - A number between -2 and 2 calculated as +/- Coverage + Synteny. Consult the manual for an explanation of this parameter.
Confidence - The confidence attached to the contig classification. This can be IDEAL / HIGH / MEDIUM / LOW / NA .  
Classification - The classification of this contig, which can either be Chromosome, the plasmid identifier followed by two digit type or Ambiguous. Consult the manual for an explanation of this parameter.
Metadata - Extra information on regarding the plasmid classification calculations. Consult the manual for an explanation of this parameter.
