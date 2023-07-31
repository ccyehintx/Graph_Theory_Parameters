# IPython log file
import matplotlib.pyplot as plt
import json
import numpy as np

f = open('input.json')
data = json.load(f)
exsnap = data['file']['extractsnap']
exincr = eval(data['file']['extract_incre'])
idx_ls = list(range(int((exsnap[1] - exsnap[0])/exincr)))

deg = np.load('deg.npy', allow_pickle=True)
bc = np.load('bc.npy', allow_pickle=True)
neg = np.load('neg.npy', allow_pickle=True)
pos = np.load('pos.npy', allow_pickle=True)
bondc = np.load('bondc.npy', allow_pickle=True)
csize = np.load('csize.npy', allow_pickle=True)

fig, ax = plt.subplots(1,1)
# Cluster size
for idx in range(len(csize)):
    MyList = csize[idx]
    my_dict = {i:MyList.count(i) for i in MyList}
    ll = list(my_dict.items())
    x = []
    y = []
    la = 'Timestep at {}'.format(idx*exincr+exsnap[0])
    for i in ll:
        x.append(i[0])
        y.append(i[1])
    #plt.scatter(x, y, label=la)
    ax.scatter(x, y, label=la, marker='.')
    
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel("Cluster Size")
ax.set_ylabel("Count")
ax.legend()
plt.show()

# Degree list
degmat = np.zeros((7, len(deg)))
for Tanals in range(len(deg)):
    anals = list(deg[Tanals])
    my_dict = {i:anals.count(i) for i in anals}
    for j in list(my_dict.items()):
        degmat[j[0]][Tanals] = j[1]
plt.clf()
for i in range(len(degmat)):
    xdeg = np.array(list(range(len(degmat[i]))))*exincr+exsnap[0]
    ydeg = degmat[i]
    ll = 'd={}'.format(i)
    plt.plot(xdeg, ydeg, label=ll)
plt.legend()
plt.xlabel('Timestep')
plt.ylabel('Counts')
plt.show()

# Bond formation/breakage
plt.clf
pp = np.round(pos)
nn = np.round(neg)
pn = np.round(pos) - np.round(neg)
xx = np.array(list(range(len(pp))))*exincr+exsnap[0]
plt.plot(xx, pp, label='New bond formation')
plt.plot(xx, nn, label='New bond breakage')
plt.plot(xx, pn, label='Total change')
plt.legend()
plt.show()
