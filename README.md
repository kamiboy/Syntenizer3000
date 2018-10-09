# Syntenizer3000
Tool for bacterial genome analysis including but not limited to synteny measurement.

Run build.sh in the folder with the source code to build an executable on a UNIX derived OS. GCC with at least C++11 support or equivalent must be installed.

Run without any input parameters to get a list of possible arguments and a short explanation of what they do.

Below is an example of running the program with a full suite of analysis enabled:

./Syntenizer --gffdir=/data/input --outdir=/data/output --genegroups=/data/input/gene_groups.csv --contigcolours=/data/input/ContigColours.csv --syntenizegroups --generatestraincharts --generatesyntenycharts --chromosomeidentifier=chromosome --plasmididentifier=Rh --fragmentidentifier=fragment --genegrouphighlightcolours=/data/input/GeneGroupHighlightColours.csv --ffndir=/data/input --generatefnagroups --presenceabsencematrix --gapitpresenceabsencematrix --generelationmatrix --classifyfragments --disambiguategroups

In the above example all input files, including *.gff, *.fnn gene groups and others are located in the folder '/data/input', and all files generated by the tool are output into the folder '/data/output'.

If the program has been configured to output circular strain chart data or gene group synteny map data then the next step is to render these using the provided R scripts 'GenerateCharts.R' and 'GenerateSyntenyMaps.R'. These scripts need to be rung with a single input parameter, which is the path to the folder with the chart or map data which in the above example would be '/data/output'.

Things to note:

If issues or uncertainties crop up while trying to use the program you can refer to the example dataset provided. A full suit of analysis can be performed using said dataset so it should serve as a guide on how input data should be formatted for the program to work.

Only GFF3 files generated by PROKKA have been tested. Other versions or files generated by other tools may contain variations that might cause incompatibility issues. In that case contact the author for possible solutions or fixes.

Two gene group file formats are supported. The direct gene group file output from ProteinOrtho, plus an internal format detailed below.
Each gene group in the internal format is defined on a ';' separated line starting with The group id as the first item. All subsequent items should be a contatanation of the strain id and gene id with the character '|' separating them. See gene_groups.csv in the example dataset for an example.

To generate basic circular charts the program only needs properly formatted GFF3 files and accompanying gene group file.
However, the GC3s lane requires the presence of gene sequence data. If these are part of the GFF3 files then they will be parsed automatically. Alternatively they can be provided using gene group *.fna files.

In order for contigs to be coloured according to their type (chromosome/plasmid/fragment) then the contig id's in the GFF3 files have to contain an string patter that identifies them as such. The program can then be configured to recognize these 

Contig colours...

To Be Completed
