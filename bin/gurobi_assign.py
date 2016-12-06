#!/usr/bin/env python

from gurobimh import *
import numpy as np
import pandas
from pandas import *
import re
import sys
from os.path import join
import scipy.sparse
from scipy.sparse import find
from scipy.sparse import coo_matrix

def sparse_df_to_array(df):
    num_rows = df.shape[0]   

    data = []
    row = []
    col = []

    for i, col_name in enumerate(df.columns):
        if isinstance(df[col_name], pandas.SparseSeries):
            column_index = df[col_name].sp_index
            #print type(column_index)
            if isinstance(column_index, pandas._sparse.BlockIndex):
                column_index = column_index.to_int_index()

            ix = column_index.indices
            data.append(df[col_name].sp_values)
            row.append(ix)
            col.append(len(df[col_name].sp_values) * [i])
        else:
            data.append(df[col_name].values)
            row.append(np.array(range(0, num_rows)))
            col.append(np.array(num_rows * [i]))

    data_f = np.concatenate(data)
    row_f = np.concatenate(row)
    col_f = np.concatenate(col)

    arr = coo_matrix((data_f, (row_f, col_f)), df.shape, dtype=np.float64)
    return arr.tocsr()


nrcols=30
col_names = range(-1,nrcols)

cur_dir=sys.argv[1]

ref_ind = pandas.read_csv(join(sys.argv[1],'genome_info.tr'),delim_whitespace=True,names=['ref','len'],skiprows=1,usecols=[1,2])
genomes = pandas.read_csv(join(sys.argv[1],'genome_info.tr'),delim_whitespace=True,header=None,nrows=1)
D1 = pandas.read_csv(join(sys.argv[1],'abundance.m_alphas'),delim_whitespace=True,names=['theta'],header=None,skiprows=7,
usecols=[0]) # skip header and noise term
data = pandas.read_csv(join(sys.argv[1],'abundance.m_phi'),delim_whitespace=True,header=None,names=col_names,index_col=0,skiprows=5).to_sparse()

numT = genomes[2][0]#np.int32(data.shape[1]/float(2))
numC = data.shape[0]

data2=sparse_df_to_array(data)

row=[]#np.empty(data2.size/2)
col=[]#np.empty(data2.size/2)
val=[]#np.empty(data2.size/2)
I,J,V=find(data2)
#for i in range(len(I)/2):

#ind=0
k=0
while k < len(I)-1:
    #if V[k]>1e-5:
    i=I[k]
    c=V[k]
    v=V[k+1]

    if not c<=0.1:
        row.append(i)
        col.append(c-1)
        val.append(v)
        
        k=k+2
    else: k=k+1

row=np.int32(row).tolist()
col=np.int32(col).tolist()

Probabilities=coo_matrix( (val, (row, col)) , shape=(numC, numT))

Probabilities=Probabilities.tocsr()
row_ind=np.split(Probabilities.indices, Probabilities.indptr[1:-1])
Probabilities=Probabilities.tocsc()
col_ind=np.split(Probabilities.indices, Probabilities.indptr[1:-1])


D1=np.array(D1['theta'])
Dist = np.int32(np.round(D1/float(np.sum(D1))*numC)).tolist()

#Dist[1]=Dist[1]-1

print np.sum(Dist)
print numC

m=Model("Assignment")
#reads=data.index.get_values()

#X = [None]*len(I)

#X = np.empty((numT,numC))
X=[x[:] for x in [[None] * numT] * numC]
for i,r_i in enumerate(row_ind):
    #X.append([None]*numC)

    #nonzero=find(col)
    for j in r_i:
        #print i,j
        X[i][j]=m.addVar(vtype=GRB.BINARY,name="x%d%d" % (i,j) )

m.update()
m.modelSense = GRB.MAXIMIZE
constraintT = []
constraintC = []
#constraintTot = []

for i,c_i in enumerate(col_ind):
    #nonzero=find(ProbabilitiesT[i])
    #nonzero=nonzero[1].tolist()
    #t=i
    #for k in nonzero:
    #    print ProbabilitiesT[i][k]
    #print nonzero
    #for p in ProbabilitiesT[i]:
    #	print p
    #[X[1][50181] for c in nonzero]
    #a0=[X[i][c] for c in nonzero[:-2]]
    #a=[X[1][50181] for c in [0,1,2]]
    #constraintT.append(m.addConstr(quicksum(a) == Dist[1],'constraintT'))
    constraintT.append(m.addConstr(quicksum(X[r][i] for r in c_i) == Dist[i],'constraintT%d' % i))

for i,r_i in enumerate(row_ind):
    constraintC.append(m.addConstr(quicksum(X[i][c] for c in r_i) == 1,'constraintC%d' % i))

#constraintTot.append(m.addConstr(quicksum(quicksum([X[r][c] for r in c_i]) for c,c_i in enumerate(col_ind)) == numC,'constraintTot'))

Probabilities=Probabilities.todok()
m.setObjective(quicksum(quicksum([X[r][c]*Probabilities[r,c] for r in c_i]) for c,c_i in enumerate(col_ind)))

m.update()
m.optimize()


out=np.zeros(numC)
matches=np.zeros(numT)
allmatch=np.zeros(numT)
reg=['']*numT

#nt ref_ind['ref']

# populate from genome info
#for i,ref in enumerate(ref_ind['ref']):
#    reg[i]=re.match('.+(.+)\.',ref).group(0)[:-1]

read_names=data.index.values

for c,c_i in enumerate(col_ind):
    for r in c_i:
        #print "     ",
        #tmp=int(X[t][c].getAttr("x"))
        if int(X[r][c].getAttr("x")):
            out[r]=c
            allmatch[c]=allmatch[c]+1
            if re.match(ref_ind['ref'][c],read_names[r]):
                matches[c]=matches[c]+1

print matches
print Dist
print np.array(Dist)/float(2)