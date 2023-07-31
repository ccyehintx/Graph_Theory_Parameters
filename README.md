# Graph_Theory_Parameters
This is to extract coordinates and compute graph theory parameters of colloidal systems
Graph theory is a useful tool to provide topological descriptions to the colloidal system by analyzing its connectivity with all other colloids
Some graph thery parameters might be useful depending on the users' area of interest

Edit the input.json file first to specify the snapshot file to read from

In the file section:

"exsnap" gives the range of time steps you wish to extract and "exincre" is the stepsize

"node_idx" and "edge_idx" are the labeling of the node(colloid) and edge(linker)

"cutoff" is the cutoff of the distance for the program to read as bonded

"patch_number" is the number of patches on a colloid

"bead_number" is the number of beads on a linker (excluding the heads)


Some useful graph theory parameters can be computed by this program by simply change to "on", otherwise "off":
"Density"
"Average_degree"
"Largest_size_of_cluster"
"Global_efficiency"
"Assortativity_coefficient"
"Average_clustering_coefficient"
"Average_nodal_coefficient"
"Betweenness_centrality"
"Closeness_cetrality"

Running the script data_adjacency.py will give some more analysis based on the use of adjacency matrix, such as cluster distribution and coordination number:
![image](https://github.com/ccyehintx/Graph_Theory_Parameters/assets/124641066/fbea2d3f-87ab-49ca-8486-3a6070a1c9c6)
![Figure_2](https://github.com/ccyehintx/Graph_Theory_Parameters/assets/124641066/28b4692d-989e-4c8b-8aab-0f63256ca621)




