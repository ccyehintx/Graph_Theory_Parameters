# This file will record the adjacency matrix for several computations in graph theory
# This will only output the files with parameters, another script is written to analyse the results (to better visualization)

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
        
###################
col_csize = []
col_deg = []
col_bc = []
col_pos = []
col_neg = []
col_bondc = []
A0 = 'empty_ls'
for nstep in range(exsnap[0], exsnap[1], exincr):
    sidx = 31009*nstep
    fidx = sidx + 31008
    sort_info = sort_info_func(mylines, sidx, fidx)
    patch_ls, head_ls = pandh_ls(sort_info, nodeidx, edgeidx)
    bond_ls = bond_ls_func(sort_info, linker_num, patch_ls, head_ls, npatch, nbead, colloid_num, bcutoff)
    realbond = []
    for i in bond_ls:
        if len(i) == 2:
            realbond.append(i)
    G = gen_G(colloid_num, realbond)
    edge_list = G.edges()
    degree_ls = nx.degree(G)
    big_cluster = list(nx.connected_components(G))
    ttcnt = []
    if len(big_cluster) == 0:
        ttcnt.append(colloid_num)
    else:
        for i in big_cluster:
            ttcnt.append(len(i))
    # The count for cluster size has completed
    dll = list(degree_ls)
    for i in range(len(dll)):
        dll[i] = dll[i][1]
    #degree in list form has cmpleted
    between_cen = nx.betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False, seed=None)
    bcls = np.array(list(between_cen.values()))
    #betweenness centrality collection completed
    Aprep = nx.adjacency_matrix(G).todense()
    A = np.squeeze(np.asarray(Aprep))
    #first step to generate A is completed
    pos_bondc = 0
    neg_bondc = 0
    ttbcall = []
    if A0 == 'empty_ls':
        A0 = A
    else:
        delA = A - A0
        A0 = A
        ttp = []
        ttpl = []
        for i in range(len(delA)):
            ttbc = 0 # this is the bond change of this colloid
            for j in range(-6, 7):
                ccc = list(np.where(delA[i] == j)[0])
                if (j % 2) == 0: # this is to remove self-loops
                    if i in ccc:
                        ccc.remove(i)
                if j > 0:
                    pos_bondc = pos_bondc + abs(j)*len(ccc)
                else: # it wont affect if j is 0 
                    neg_bondc = neg_bondc + abs(j)*len(ccc)
                ttbc = ttbc + abs(j)*len(ccc)
            ttbcall.append(ttbc)
    col_csize.append(ttcnt)
    col_deg.append(dll)
    col_bc.append(bcls)
    col_pos.append(pos_bondc/2) #divided by 2 as overcounting
    col_neg.append(neg_bondc/2)
    col_bondc.append(ttbcall)
    print('Info at {} is collected'.format(nstep))

    
np.save('csize.npy', np.array(col_csize))
np.save('deg.npy', np.array(col_deg))
np.save('bc.npy', np.array(col_bc))
np.save('pos.npy', np.array(col_pos))
np.save('neg.npy', np.array(col_neg))
np.save('bondc.npy', np.array(col_bondc))
# ttpl is the total count for the plus 1 bond, need to do for other numbers later
# here can also try to remove the diagonal term which is the self loop
# then add up all the counts so we can get the total bond change count
