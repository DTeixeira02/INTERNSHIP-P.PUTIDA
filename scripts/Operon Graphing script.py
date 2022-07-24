#Graphing_script-input:IGR_Operon_assoc
import networkx as nx
import os
import re
from collections import defaultdict
import csv
import random
import matplotlib.pyplot as plt

def parse_igr_pred(igrdoc, genid):
    igr_data = []
    igr_dict = defaultdict(list)
    with open(igrdoc, "r") as file:
        for line in file.readlines():
            currentline = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","")  
            currentline = tuple(currentline.split(","))  
            if currentline[1] == genid and len(currentline) == 3:  
                igr_data.append(currentline)
            elif len(currentline) != 3:  
                print(f"Potential data issue, {line=}")
    for igr, genome, gene in igr_data:
        igr_dict[igr].append(gene)
    return igr_dict
            
def parse_op_pred(opdoc):
    op_data = []
    with open(opdoc, "r") as file1:
        for line in file1.readlines():
            currentline = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","")  
            currentline = tuple(currentline.split(","))
            if len(currentline) == 2:
                op_data.append(currentline)
            else:
                print(f"Potential data issue, {line=}")
    operon_dict = {gene: opn for opn, gene in op_data}
    return operon_dict

datasetprompt = input("Which dataset? (comp/both):")
if datasetprompt == "comp":
    igrloc = "GENOMES\\piggy_out_BOTHGEN\\IGR-Gene associations per genome(BOTHGENOMES).txt"
    oploc = "GENOMES\\operon_predictions_comp"
    graphloc = "GENOMES\\GraphML_compgen\\IGR_OPERON_GML-"
elif datasetprompt == "both":
    igrloc = "GENOMES\\piggy_out_BOTHGEN\\IGR-Gene associations per genome(BOTHGENOMES).txt"
    oploc = "GENOMES\\operon_predictions_both"
    graphloc = "GENOMES\\GraphML_bothgen\\IGR_OPERON_GML-"

collection = os.listdir(oploc)

for i in range(len(collection)):
    currentgenomeloc = str(oploc + "\\" + collection[i])
    currentgenomeid = str(re.findall("Operon_prediction_(.*).txt", collection[i])[0])
    igrdict = parse_igr_pred(igrloc, currentgenomeid)
    opdict = parse_op_pred(currentgenomeloc)
    igr_op = defaultdict(list)
    for igr, genes in igrdict.items():
        for gene in genes:
            try:
                igr_op[igr].append(opdict[gene])  
            except KeyError:
                if gene != "N/A-IGRin/isdoubleterminator":  
                    print(f"IGR {igr} is not a double-terminator, but no operon found for {gene}")
    igr_graph = nx.DiGraph()
    for igr, operons in igr_op.items():
        for op in operons:
            igr_graph.add_edge(igr, op)
    nx.write_graphml(igr_graph, str(graphloc + currentgenomeid))  
    print(f"The task is {(i/len(collection))*100} percent done")
