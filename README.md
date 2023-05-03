# rg3m
##### Adam Podgorny & Cory Jenkinson
###### Care of: Oakley Lab, University of Kansas
<br>


### Overview
This script takes in a query and target file, and then returns matches within a specified cutoff. 

### Algorithm

The input consists of a BLAST formatted resistance gene file, secondary metabolite file, and a cutoff distance in base-pairs. The optional 'nocheck' flag affects filtering behavior for the query genes. Additionally, a homolog search mode is available for finding resistance gene homologs for successful searches.

The resistance gene file is read in, such that a table is created. Each entry contains the gene name, chromosome, chromosomal start and end locations, E-Value, organism name and organism ID. The center of the gene is calculated by the simple geometric average of the start and end locations. If the 'nocheck' parameter is not passed, genes names without duplicates are filtered out. 

The secondary metabolite file is then read in to create a gene table as above (without any filtering by redundancy).

Each entry in the resistance gene table is then tested against each entry in the secondary metabolite gene table. If the organism names match, the entries are checked for a chromsome match. If the chromosome matches, genomic overlap and inclusion within the cutoff are then assessed. If the putative hit is within the cutoff, and the two genes do not overlap, the hit is accepted and output.

The accepted hits are then printed to screen. If an output file is specified, these entries are output in csv format. Each hit entry contains the organism name, organism ID, chromosome, resistance gene start, resistance gene end, resistance gene E-value, secondary metabolite chromosome (SM), secondary metabolite gene start, secondary metabolite gene end, and the distance between genes.

If the homolog mode is activated, any resistance genes that match to a hit will be used to rescan the resistance gene file for homologs in the same organism. This will be output to the specified filename.


<br>


### Usage

`python functional_gene_intersection --cutoff <cutoff in bp> --resistance_gene <resistance gene file> --sm_gene <SM gene file> [--out output csv] [--rghomologs resistance gene homolog mode] [--nocheck]`
<br>

The Python version is up to choice, although this was developed on a python2 base. It requires no libraries beyond the libraries that come standard with Python2.x and 3.x to make deployment simpler.

* `--cutoff` - This parameter controls the distance from the query 
cluster center from which hits are accepted or rejected. The units are nts.

* `--resistance_gene` - The file with the resistance gene files.

* `--sm_gene` -- The file with the secondary metabolite genes.

* `--out` - An optional parameter which will allow the user to dump the results to the specified csv.

* `--nocheck` - An optional parameter which disables exclusion of duplicate gene names from the resistance genes.

* `--rghomologs` - An optional parameter which turns on homolog mode for the resistance genes.

* `--gene_length_cutoff` - An optional parameter which controls the maximum gene size. Genes above the threshold are excluded. Default 50,000 nts.
