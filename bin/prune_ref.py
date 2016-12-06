#!/usr/bin/env python
import sys
import newick
import numpy as np

def map_tree(obj,fun,*args):
    if type(obj)==list:
        fun(obj,*args)
        return [map_tree(node,fun,*args) for node in obj]
    elif type(obj)==newick.Node:
        if not obj.descendants:
            return fun(obj,*args)
        else:
            return map_tree(obj.descendants,fun,*args)

def prune(obj,threshold):
    if type(obj)==list:
        for sub_obj in obj:
            if len(sub_obj.descendants)==1 and sub_obj.length:
                sub_obj.length=sub_obj.length+sub_obj.descendants[0].length
                sub_obj.descendants=[]
        else:
            removed=True
            while(removed):
                removed=False
                # dont touch nodes with descendants
                touch_ind=[]
                for i in range(len(obj)):
                    if not obj[i].descendants: 
                        touch_ind.append(i)

                # sort nodes according to branch length
                lengths=[node.length for node in obj]
                pair_ind=np.triu_indices(len(touch_ind), 1)
                
                for k in range(len(pair_ind[0])):
                    i,j=pair_ind[0][k],pair_ind[1][k]
                    ll=[lengths[touch_ind[i]],lengths[touch_ind[j]]]
                    tot=sum(ll)
                    if tot<threshold:
                        ind=ll.index(min(ll))
                        removed=True
                        if ind==0: del obj[touch_ind[i]]
                        else: del obj[touch_ind[j]]
                        break

def get_names(obj):
    if type(obj)==newick.Node: print obj.name

def get_lengths(obj):
    if type(obj)==newick.Node: print obj.length

if __name__ == "__main__":
	threshold = 0.001
	tree = newick.read(sys.argv[1])

	map_tree(tree,prune,threshold)
	map_tree(tree,prune,threshold)
	map_tree(tree,get_names)

	newick.write(tree,sys.argv[1]+'_pruned')