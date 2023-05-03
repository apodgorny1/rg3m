# **FUN**ctional **G**ene **I**ntersection
##### Adam Podgorny & Cory Jenkinson
###### Care of: Oakley Lab, University of Kansas
<br>


### Overview
This script takes in a query and target file, and then returns matches within a specified cutoff. 

### Algorithm

The input consists of a BLAST formatted query gene file, target gene file, and a cutoff distance in base-pairs. The optional 'nocheck' flag affects filtering behavior for the query genes.

The query file is read in, such that a table is created. Each entry contains the gene name, chromosome, chromosomal start and end locations, E-Value, organism name and organism ID. The center of the gene is calculated by the simple geometric average of the start and end locations. If the 'nocheck' parameter is not passed, genes names without duplicates are filtered out. 

The target file is then read in to create a gene table as above (without any filtering by redundancy).

Each entry in the query gene table is then tested against each entry in the target gene table. If the organism names match, the entries are checked for a chromsome match. If the chromosome matches, genomic overlap and inclusion within the cutoff are then assessed. If the putative hit is within the cutoff, and the two genes do not overlap, the hit is accepted and output.

The accepted hits are then printed to screen. If an output file is specified, these entries are output in csv format. Each hit entry contains the organism name, organism ID, chromosome, resistance gene start, resistance gene end, resistance gene E-value, secondary metabolite chromosome (SM), secondary metabolite gene start, secondary metabolite gene end, and the distance between genes.


<br>


### Usage

`python functional_gene_intersection _cutoff_ _query_file.csv_ _target_file.csv_ _[outputfilename.csv - optional]_ _[nocheck - optional]_`  
<br>

The Python version is up to choice, although this was developed on a python2 base.

* `cutoff` - This parameter controls the distance from the query 
cluster center from which hits are accepted or rejected. The units are bp

* `query_file.csv` - The filename with the query clusters

* `target_file.csv` - The filename with the target genes

* `outputfilename.csv` - An optional parameter which will allow the user to dump the results to a csv directly

* `nocheck` - An optional parameter which disables exclusion of duplicat gene names from the hit result


