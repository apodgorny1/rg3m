#Adam Podgorny & Cory Jenkinson
import os
import argparse

class gene:
    """Storage class to hold resistance gene or secondary metabolite gene descriptions."""
    def __init__(self, start_, end_, scaffold_, org_name_ = "", evalue_ = -1, pct_id_ = -1):
        self.start = float(start_)
        self.end = float(end_)
        self.scaffold = scaffold_
        self.evalue = evalue_
        self.org_name = org_name_
        self.pct_id = pct_id_
        self.average = (self.start + self.end + 0.0)/2.0 ##The coordinates of a gene are defined from the mean.
    def queryout(self):
        """Helper class to output vital resistance gene statistics."""
        if (self.pct_id == -1):
            output = str(self.org_name + ", " + self.scaffold + ", " + str(int(self.start)) + ", " + str(int(self.end)) + ", " + str(self.evalue))
        else:
            output = str(self.org_name + ", " + self.scaffold + ", " + str(int(self.start)) + ", " + str(int(self.end)) + ", " + str(self.evalue) + ", " + str(self.pct_id))
        return output
    def targetout(self):
        """Helper class to output vital secondary metabolite statistics."""
        output = str(str(self.scaffold) + ", " + str(int(self.start)) + ", " + str(int(self.end)))
        return output


parser = argparse.ArgumentParser(description="Checks a query file line by line against a query database.")
parser.add_argument("--cutoff", required=True, help = "Distance cutoff for eligible targets in nucleotides.")
parser.add_argument("--resistance_gene", required=True, help = "Resistance gene database file.")
parser.add_argument("--sm_gene", required=True, help = "Secondary metabolite database file.")
parser.add_argument("--out", help="Optional: File into which to write output.")
parser.add_argument("--nocheck", action='store_true', help="If flag is active, skips the duplicate check.")
parser.add_argument("--gene_length_cutoff", help="Gene Length Cutoff. Default = 50,000 bp", default=50000)
parser.add_argument("--rghomologs", help="Option: Print out the homologous genes in the resistance gene db? If so, to where")


###Argument Handling
args = parser.parse_args()

cutoff = int(args.cutoff)
query_file = args.resistance_gene
target_file = args.sm_gene
GENE_LENGTH_CUTOFF = int(args.gene_length_cutoff)

if (args.rghomologs):
    h_flag = True
    h_out = args.rghomologs
else:
    h_flag = False 


if (args.nocheck):
    nc_flag = True
else:
    nc_flag = False    

if (args.out):
    outflag = True
    outfile = args.out
else:
    outflag = False
    
##NOTE: The target file needs the Location column removed or else it breaks the split

###Opening the files
query_fd = open(query_file)
target_fd = open(target_file)



queries = []
targets = []
print("\nUsing cutoff of " + str(cutoff) + " bp.")

##Read in the resistance gene file

##strip the header out
line_in = query_fd.readline()
line_in = line_in.split(",")
t_line_in = []
for x in line_in:
    t_line_in.append(x.strip())
    
line_in = t_line_in

query_header_len = len(line_in)
if (line_in[0][0] == "\""):
 l_t = [x[1:-1] for x in line_in]
 line_in = l_t
if ("Hit Start") in line_in: 
    query_start_c = line_in.index("Hit Start")
elif ("hit start") in line_in: 
    query_start_c = line_in.index("hit start")
elif ("hit_start") in line_in: 
    query_start_c = line_in.index("hit_start")
elif ("Start") in line_in: 
    query_start_c = line_in.index("Start")
elif ("start") in line_in:
    query_start_c = line_in.index("start")

#We know the end will always be one after
query_end_c = 1 + query_start_c

if "EValue" in line_in:
    e_value_c = line_in.index("EValue")
elif "evalue" in line_in:
    e_value_c = line_in.index("evalue")
elif "Evalue" in line_in:
    e_value_c = line_in.index("Evalue")
else:
    e_value_c = -1
    
#Organism name
if "Organism" in line_in:
    org_name_c = line_in.index("Organism")
elif "organism" in line_in:
    org_name_c = line_in.index("organism")
elif "Organism Name" in line_in:
    org_name_c = line_in.index("Organism Name")
elif "\"Organism Name\"" in line_in:
    org_name_c = line_in.index("\"Organism Name\"")
elif "\"Organism\"" in line_in:
    org_name_c = line_in.index("\"Organism\"")
else:
    print("Malformed Resistance Gene Header. Please rename the organism name column to 'Organism'\n")
    exit()


if "scaffold" in line_in:
    scaffold_c = line_in.index("scaffold")
elif ("Scaffold") in line_in: 
    scaffold_c = line_in.index("Scaffold")
elif "Chromosome" in line_in:
    scaffold_c = line_in.index("Chromosome")
elif "chromosome" in line_in:
    scaffold_c = line_in.index("chromosome")
elif ("Hit Name") in line_in:
    scaffold_c = line_in.index("Hit Name")
else:
    scaffold_c = -1

if "Location" in line_in:
    location_c = line_in.index("Location")
elif ("location") in line_in: 
    location_c = line_in.index("location")
else:
    location_c = -1

if "% Identity" in line_in:
    id_pct_c = line_in.index("% Identity")
elif ("% identity") in line_in: 
    id_pct_c = line_in.index("% identity")
elif ("Percent Identity") in line_in:
    id_pct_c = line_in.index("Percent Identity")
elif "Identity Percent" in line_in:
    id_pct_c = line_in.index("Identity Percent")
elif "Per. Identity" in line_in:
    id_pct_c = line_in.index("Per. Identity")
elif "Per Identity" in line_in:
    id_pct_c = line_in.index("Per Identity")    
elif "Per Ident" in line_in:
    id_pct_c = line_in.index("Per Ident")    
elif "Identity" in line_in:
    id_pct_c = line_in.index("Identity")     
elif "Identities" in line_in:
    id_pct_c = line_in.index("Identities")
elif "% Hit Identity" in line_in:
    id_pct_c = line_in.index("% Hit Identity")      
else:
    id_pct_c = -1

###Iterate over the resistance gene entries
for i in query_fd:
  ###Grab each line and split it into the necessary fields.
  ###There is some text handling due to some observed cases in each field.
  line_in = i
  if (location_c != -1):
   if (" (+)" in line_in):
    cut_off = line_in.index("(+)")
   elif (" (-)" in line_in):
    cut_off = line_in.index("(-)")
   else:
    print("Nonstandard location format. Please delete the location column in the resistance gene db file\n")
    exit(1)
   cut_off = cut_off + 3
   line_in = line_in[0:location_c] + line_in[cut_off:len(line_in)]
  row = line_in.split(",")
  hit_start = row[query_start_c]
  hit_end = row[query_end_c]
  scaffold = row[scaffold_c].strip()
  if (scaffold[0] == "\""):
   scaffold = scaffold[1:len(scaffold) - 1].strip()
  org_name = row[org_name_c]
  org_name = org_name.replace("\"", '').strip()


  if (e_value_c != -1):
   e_value = row[e_value_c]
   e_value = e_value[1:len(e_value)]
  else:
    ##Handle cases without evalues
    e_value = 0.0

  if (id_pct_c != -1):
   pct = row[id_pct_c]
  else:
   pct = -1
   


  if (len(row) != query_header_len):
       ###Handle cases when a field is omitted for an entry
	   q_row_fallback_counter = 0
	   for j in row:
		   if (j == ("(+)") or (j == "(-)")):
			   hit_end = row[(q_row_fallback_counter-1)]
			   hit_start = row[(q_row_fallback_counter - 2)]
			   scaffold = row[(q_row_fallback_counter - 3)]
			   if (scaffold[0] == "\""):
				   scaffold = scaffold[1:len(scaffold) - 1]
		   q_row_fallback_counter = q_row_fallback_counter + 1

  queries.append(gene(hit_start, hit_end, scaffold, org_name, evalue_=e_value, pct_id_=pct))



print("Resistance Gene File Imported")

if (nc_flag == False):
 ##Fixing duplicates in said query file
 s_t = [x.org_name for x in queries]
 ##Checking Gene length here now too
 q_temp = [x for x in queries if s_t.count(x.org_name) > 1 if (abs(x.end - x.start) <= GENE_LENGTH_CUTOFF)]
 queries_old = queries
 queries = q_temp


##Now for the Secondary Metabolites
line_in = target_fd.readline()
line_in = line_in.split(",")
t_line_in = []
for x in line_in:
    t_line_in.append(x.strip())
    
line_in = t_line_in
target_header_len = len(line_in)
if (line_in[0][0] == "\""):
 l_t = [x[1:-1] for x in line_in]
 line_in = l_t

if ("Hit Start") in line_in: 
    target_start_c = line_in.index("Hit Start")
elif ("hit start") in line_in: 
    target_start_c = line_in.index("hit start")
elif ("hit_start") in line_in: 
    target_start_c = line_in.index("hit_start")
elif ("Start") in line_in: 
    target_start_c = line_in.index("Start")
#we know the end will always be one after
target_end_c = 1 + target_start_c

if "EValue" in line_in:
    t_e_value_c = line_in.index("EValue")
elif "evalue" in line_in:
    t_e_value_c = line_in.index("evalue")
elif "Evalue" in line_in:
    t_e_value_c = line_in.index("Evalue")
else:
    t_e_value_c = -1
    
#organism name
if "Organism" in line_in:
    t_org_name_c = line_in.index("Organism")
elif "organism" in line_in:
    t_org_name_c = line_in.index("organism")
elif "Organism Name" in line_in:
    t_org_name_c = line_in.index("Organism Name")
elif "\"Organism Name\"" in line_in:
    t_org_name_c = line_in.index("\"Organism Name\"")
elif "\"Organism\"" in line_in:
    org_name_c = line_in.index("\"Organism\"")
else:
    print("Malformed Resistance Gene Header. Please rename the organism name column to 'Organism'\n")
    exit()


if "scaffold" in line_in:
    t_scaffold_c = line_in.index("scaffold")
elif ("Scaffold") in line_in: 
    t_scaffold_c = line_in.index("Scaffold")
elif ("Hit Name") in line_in:
    t_scaffold_c = line_in.index("Hit Name")
elif "Chromosome" in line_in:
    t_scaffold_c = line_in.index("Chromosome")
elif "chromosome" in line_in:
    t_scaffold_c = line_in.index("chromosome")
else:
    t_scaffold_c = -1

if "Location" in line_in:
    location_c = line_in.index("Location")
elif ("location") in line_in: 
    location_c = line_in.index("location")
else:
    location_c = -1
    
    
##Strip headers again
bc = 0
for i in target_fd:
  bc = bc + 1
  line_in = i
  if (location_c != -1):
   if ("(+)" in line_in):
    cut_off = line_in.index("(+)")
   elif ("(-)" in line_in):
    cut_off = line_in.index("(-)")
   else:
    print("Nonstandard location format. Please delete the location column in the secondary metabolite db file\n")
    exit(1)

   temp_c = line_in[cut_off]
   temp_p = cut_off
   while (temp_c != '_'):
	   temp_p = temp_p - 1
	   temp_c = line_in[temp_p]
   while (temp_c != ','):
	   temp_p = temp_p - 1
	   temp_c = line_in[temp_p]
   cut_off = cut_off + 1
   line_in = line_in[0:temp_p] + line_in[cut_off:len(line_in)]
  row = line_in.split(",")
  org_name = row[t_org_name_c]
  org_name = org_name.replace("\"", '').strip()
  scaffold = row[t_scaffold_c].strip()
  if (scaffold[0] == "\""):
   scaffold = scaffold[1:len(scaffold) - 1].strip()
  hit_start = row[target_start_c]
  hit_end = row[target_end_c]

  if (len(row) != target_header_len):
	   q_row_fallback_counter = 0
	   for j in row:
		   if (j == ("(+)") or (j == "(-)")):
			   hit_end = row[(q_row_fallback_counter-1)]
			   hit_start = row[(q_row_fallback_counter - 2)]
			   scaffold = row[(q_row_fallback_counter - 3)]
			   if (scaffold[0] == "\""):
				   scaffold = scaffold[1:len(scaffold) - 1]
		   q_row_fallback_counter = q_row_fallback_counter + 1


  targets.append(gene(hit_start, hit_end, scaffold, org_name, ""))


target_fd.close()
query_fd.close()

print("Secondary Metabolite File Imported")
print("Accepted Resistance Genes: " + str(len(queries)) + "\n")
q_length = len(queries)
t_length = len(targets)
print(str(q_length) + " Resistance Genes")
print(str(t_length) + " Secondary Metabolites\n")
print("Results:")
print("__________________")


##Header for writing output to a csv.
if (id_pct_c == -1):
    outheader = "Organism Name, Resistance Gene Scaffold, Resistance Gene Start, Resistance Gene End, Resistance Gene E-value, SM Gene Scaffold, SM Gene Start, SM Gene End, Distance\n"
else: 
    outheader = "Organism Name, Resistance Gene Scaffold, Resistance Gene Start, Resistance Gene End, Resistance Gene E-value, Resistance Gene % Identity, SM Gene Scaffold, SM Gene Start, SM Gene End, Distance\n"

if (outflag == True):
    ##Create an output file if flag is active
    fdout = open(outfile, "w")
    fdout.write(outheader)

if (h_flag == True):
    h_confirms =[]
    h_fd = open(h_out, 'w')
    if (id_pct_c != -1):
        h_fd.write("Organism Name, E-value, % Identity, Homologous Gene Scaffold, Homologous Gene Start, Homologous Gene End\n")
    else:
        h_fd.write("Organism Name, E-value, Homologous Gene Scaffold, Homologous Gene Start, Homologous Gene End\n")

print(outheader)

    
for i in range(q_length): ##Iterate over queries
    for j in range(t_length): ##Iterate over targets
     if (queries[i].org_name == targets[j].org_name): ##If organisms match, find clusters
      if (queries[i].scaffold == targets[j].scaffold): ##Make sure they're on same chromosome
           distance = (abs(queries[i].average - targets[j].average))
           distance_check = (distance <= cutoff) ##check distances
           ##Overlapping sequences check
           a = (queries[i].start >= targets[j].start) and (queries[i].start <= targets[j].end)
           b = (queries[i].end >= targets[j].start) and (queries[i].start <= targets[j].end)
           overlap = a and b ##Check if the genes overlap
           gene_length_check = abs(targets[j].end - targets[j].start) <= GENE_LENGTH_CUTOFF
           ##Check that the target gene is of a reasonable size
           if (not(overlap) and distance_check and gene_length_check):
            q = queries[i].queryout()
            t = targets[j].targetout()
            ##If all the conditions are met, we declare a hit, and print/write it out.
            if (outflag == True):
              nt = t + ", " + str(int(distance)) + "\n"
              fdout.write(str(q) + ", " + str(nt))
            print(str(q) + ", " + str(t) + ", " + str(int(distance)))
            
            ##If homolog mode is active, append the relevant data to the storage list
            if (h_flag == True):
                if (id_pct_c != -1):
                    homolog_list = [(h.org_name, h.evalue, h.pct_id, h.scaffold, h.start, h.end) for h in queries if (h.org_name == targets[j].org_name)] ##Find unique matches to the query
                else:
                    homolog_list = [(h.org_name, h.evalue, h.scaffold, h.start, h.end) for h in queries if (h.org_name == targets[j].org_name)] ##Find unique matches to the query
                h_confirms = h_confirms + homolog_list ##Append these to the main homolog set

##If we want to track homologs, output all of the unique resistance gene homologs.
if (h_flag == True):
    h_confirms = list(set(h_confirms)) ##
    for h in h_confirms:
        if (id_pct_c != -1):
            h_string = "{},{},{},{},{},{}\n".format(h[0], h[1], h[2], h[3], int(h[4]), int(h[5]))
        else:
            h_string = "{},{},{},{},{}\n".format(h[0], h[1], h[2], int(h[3]), int(h[4]))
        h_fd.write(h_string)


        
#Bookkeeping

if (outflag == True):
    fdout.close()
   
if (h_flag == True):
    h_fd.close()

print("\n")
