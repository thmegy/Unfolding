#!/usr/bin/python3.7

import numpy as np
import matplotlib.pyplot as plt

ngen = 25000
truth = np.random.normal(40,5,ngen)

#truth = np.append(truth, np.random.uniform(0,120,ngen*4))
truth = np.append(truth, np.random.normal(75,10,ngen//2))

#reco = np.random.normal(truth, 10, ngen)

reco = np.random.normal(truth, 0.15*truth)

#for i in truth:
#    reco = np.append(reco,np.random.normal(i, 0.15*i))

nbin = 50
bin_edges = np.linspace(0,120,nbin+1)

plt.hist(truth, bins=bin_edges, density=False, histtype='step', label='truth', color='red')
plt.xlabel('some observable')
plt.ylabel('Events/{}'.format(bin_edges[nbin]/nbin))
plt.legend()
plt.savefig('truth.pdf')
plt.close()

(bincontent, recobins, patches) = plt.hist(reco, bins=bin_edges, density=False, histtype='step', label='reco', color='green')
plt.xlabel('some observable')
plt.ylabel('Events/{}'.format(bin_edges[nbin]/nbin))
plt.legend()
plt.savefig('reco.pdf')
plt.close()

databins = []
for i in range(bin_edges.size-1):
    databins.append((bin_edges[i+1]+bin_edges[i])/2)

data = np.random.poisson(bincontent)
    
plt.hist(databins, bins=bin_edges, weights=data, density=False, histtype='step', label='data', color='black')    

plt.xlabel('some observable')
plt.ylabel('Events/{}'.format(bin_edges[nbin]/nbin))
plt.legend()
plt.savefig('data.pdf')
