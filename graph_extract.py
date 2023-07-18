# Faster way to calculate coordinates by matrix operations

import numpy as np
import gzip
import matplotlib.pyplot as plt
import random
import networkx as nx
import json

def dis_mat(coor_patch, coor_heads):
    A = np.ones((len(coor_heads), 3))
    A[:, 0] = A[:, 0]*coor_patch[0]
    A[:, 1] = A[:, 1]*coor_patch[1]
    A[:, 2] = A[:, 2]*coor_patch[2]
    B = coor_heads
    C = B - A
    D = np.square(C)
    E = np.sum(D, axis=1)
    F = np.sqrt(E)
    return F

def raw_line(line):
    info = []
    info.append(eval(line.split()[0]))
    info.append(eval(line.split()[2]))
    coor = np.array([eval(i) for i in line.split()[3:6]])
    info.append(coor)
    return info

def sort_info_func(mylines, sidx, fidx):
    sort_info = []
    for i in range(sidx+9, fidx+1):
        sort_info.append(raw_line(mylines[i]))
    return sort_info

def pandh_ls(sort_info, nodeidx, edgeidx):
    patch_ls = []
    head_ls = []
    for i in range(len(sort_info)):
        if sort_info[i][1] == nodeidx: # 2 is patch
            patch_ls.append(i) 
        elif sort_info[i][1] == edgeidx: # 3 is head
            head_ls.append(i)
    return patch_ls, head_ls

def bond_ls_func(sort_info, linker_num, patch_ls, head_ls, npatch, nbead, colloid_num, bcutoff):
    bond_ls = []
    tscol = npatch + 1
    ttcol = (npatch + 1)*colloid_num
    ttlnk = nbead + 2
    for b in range(linker_num):
        bond_ls.append(list())
    head_all_info = []
    head_iidx = []
    for j in head_ls:
        head_all_info.append(sort_info[j][2])
        head_iidx.append(sort_info[j][0])
    for i in patch_ls:
        cidx = (sort_info[i][0] - 1)//tscol
        icoor = sort_info[i][2]
        collect_dis = dis_mat(icoor, head_all_info)
        test_list = list(collect_dis)
        res = [idx for idx, val in enumerate(test_list) if val < bcutoff]
        if len(res) > 0:
            for k in res:
                bbidx = head_iidx[k]
                jjidx = (bbidx - ttcol - 1)//ttlnk #(1000 + 6*1000 = 7000)
                bond_ls[jjidx].append(cidx)
    return bond_ls


def gen_G(size, bonding):
    # Generate a graph G
    G = nx.Graph()
    for i in range(size):
        G.add_node(i)
    for i in bonding:
        G.add_edge(i[0],i[1])
    return G

def avg_deg(degree_ls):
    # Generate average degree
    sum_deg = []
    for i in degree_ls:
        sum_deg.append(i[1])
    return sum(sum_deg)/len(sum_deg)


# Read .json file (the file part)
f = open('input.json')
data = json.load(f)
gzfile = data['file']['filename']
exsnap = data['file']['extractsnap']
exincr = eval(data['file']['extract_incre'])
nodeidx = eval(data['file']['node_idx'])
edgeidx = eval(data['file']['edge_idx'])
bcutoff = eval(data['file']['cutoff'])
linker_num = eval(data['file']['linker_number'])
colloid_num = eval(data['file']['colloid_number'])
npatch = eval(data['file']['patch_number'])
nbead = eval(data['file']['bead_number'])

mylines = []
with gzip.open(gzfile, mode ='r')as file:
    for ll in file:
        mylines.append(ll)

########### Main code #################
# Read .json file (the GT paramters part)
gt_cont = []
gt_requ = []
gt_cont.append('Timestep')
for i in data['GT_para']:
    if data['GT_para'][i] == 'on':
        gt_cont.append(i)
        gt_requ.append(1)
    else:
        gt_requ.append(0)

print(gt_cont)
total_para = []
total_ncon = []
total_betc = []
total_cloc = []

for nstep in range(exsnap[0], exsnap[1], exincr):
    this_gt = []
    this_gt.append(nstep)
    sidx = 31009*nstep
    fidx = sidx + 31008
    sort_info = sort_info_func(mylines, sidx, fidx)
    patch_ls, head_ls = pandh_ls(sort_info, nodeidx, edgeidx)
    bond_ls = bond_ls_func(sort_info, linker_num, patch_ls, head_ls, npatch, nbead, colloid_num, bcutoff)
    # And this bond_ls tells all the bonding information
    # Now we have to create with linkers with both sides bonded
    realbond = []
    for i in bond_ls:
        if len(i) == 2:
            realbond.append(i)
    G = gen_G(colloid_num, realbond) # Generate a graph
    edge_list = G.edges() # List of all edges
    degree_ls = nx.degree(G) # List of degrees of all nodes
    big_cluster = list(nx.connected_components(G)) #All the clusters, if return nothing then all nodes are connected
    bigc_len = []
    for i in big_cluster:
        bigc_len.append(len(i))
    cc = nx.clustering(G)
    sum_cc = []
    for i in range(len(cc)):
        sum_cc.append(cc[i])
    p = nx.shortest_path(G) #this compute the shortest path
    if gt_requ[0] == 1: # Density
        this_gt.append(2*len(edge_list)/colloid_num/(colloid_num-1))
    if gt_requ[1] == 1: # Average degree
        this_gt.append(avg_deg(degree_ls))
    if gt_requ[2] == 1: # Largest size of cluster
        this_gt.append(max(bigc_len))
    if gt_requ[3] == 1: # Global efficiency
        this_gt.append(nx.global_efficiency(G))
    if gt_requ[4] == 1: # Assortativity coefficient
        this_gt.append(nx.degree_assortativity_coefficient(G))
    if gt_requ[5] == 1: # Average clustering coefficent
        this_gt.append(sum(sum_cc)/len(sum_cc))
    if gt_requ[6] == 1: # Nodal connectivity
        #this_gt.append(nx.node_connectivity(G))
        nconn = nx.node_connectivity(G)
        total_ncon.append(np.array(nconn))
    if gt_requ[7] == 1: # Betweenness centrality
        between_cen = nx.betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False, seed=None)
        total_betc.append(np.array(between_cen))
    if gt_requ[8] == 1: # Closeness centrality
        close_cen = nx.closeness_centrality(G)
        total_cloc.append(np.array(close_cen))
    total_para.append(np.array(this_gt))

# Storing them into .npy files
if len(total_para) > 0:
    np.save('GT_parameters.npy', np.array(total_para))
if len(total_ncon) > 0:
    np.save('GT_nodal_connectivity.npy', np.array(total_ncon))
if len(total_betc) > 0:
    np.save('GT_betweenness_centrality.npy', np.array(total_betc))
if len(total_cloc) > 0:
    np.save('GT_closeness_centrality.npy', np.array(total_cloc))
