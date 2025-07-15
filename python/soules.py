"""
=================================================================================================
 
  Copyright (C) 2025 Francois G. Meyer <FMeyer@Colorado.Edu>
 
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  
    a. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
  
    b. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
  
    c. Neither the name of the copyright holders nor the names of any
       contributors to this software may be used to endorse or promote products
       derived from this software without specific prior written permission.
  
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS on an
  "AS IS" basis. THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO
  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  BY WAY OF EXAMPLE, BUT
  NOT LIMITATION, THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO AND
  DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR
  ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE
  ANY THIRD PARTY RIGHTS.
  
  THE COPYRIGHT HOLDERS AND CONTRIBUTORS SHALL NOT BE LIABLE TO LICENSEE OR
  ANY OTHER USERS OF THIS SOFTWARE FOR ANY INCIDENTAL, SPECIAL, OR
  CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR
  OTHERWISE USING THIS SOFTWARE, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED
  OF THE POSSIBILITY THEREOF.
"""
import numpy as np
from generate_graph import genSBM
import util_funcs
np.random.seed(17)

def build_vector(x, i0, i1, k):
    """
    _Module_Name : buildvector.m
    
    _Description : construct a Soules vector. 
                  
    input: x vector of length sizex

           [i0,i1] set of indices of the block from which we construct the basis vector
           
           k: index where we split: [i0,k] on one side and [k+1,i1] on the other side

    output: the basis vector corresponding to the node


    _References :
    
    _Remarks : None
    
    _Author :                 Francois G. Meyer
    
    _Revisions History: 2025 Initial keying
    """
    # error check magnitude of k
    if (k< i0) or (k>= i1):
        print(f"\nerror: index of split {k} is out of range: {i0} to {i1}.\n")
        return None

    # potential error site w indexing
    node = x[i0-1:i1]
    #           +--------+  
    #    left:  |        |
    #           |        +---------------------
    #           i0       k                   i1
    #                    +--------------------+
    #    right:          |                    |
    #           ---------+          
    left = np.concatenate((x[i0-1:k], np.zeros(i1 - k)))
    right = np.concatenate((np.zeros(k-i0+1), x[k:i1]))
    nrm = np.linalg.norm(node)
    nrm_left = np.linalg.norm(left)
    nrm_right = np.linalg.norm(right)
    v = np.zeros(np.shape(x))
    #           +--------+  
    #           |        |
    #     v: __ i0 _____ k __________________ i1 ___
    #                    |                    |
    #                    +--------------------+
    v[i0-1:i1] = (1/nrm) * ((nrm_right/nrm_left)*left - (nrm_left/nrm_right)*right)

    return v



def topdown(M, r1):
    """
    _Module_Name : topdown
    
    _Description : explores the binary tree of Soules bases, from the top (coarsest)
                   level, down to the finest level. The function computes the best Soules
                   basis , which provides the sparest expansion of the matrix A.
    
                   The core function is buildvector that constructs a new Soules vector, based
                   on entries in r1, and the geometry in the Soules tree that is being split.
    
                   Input:
                           M the sample mean adjacency matrix, of size nNodes x nNodes
                           r1 first Soules vector: nonnegative vector; typically r1 = n^{-1/2} * one
                   Output:
                           Qn  orthogonal matrix to store all the basis vectors 
                               only nNodes - 1 columns are non zero in Qn; the 
    
    _References :
    
    _Remarks : None
    
    _Author :                 Francois G. Meyer
    
    _Revisions History: 2025 Initial keying
    
    ==============================================================================
    """
    # adjacency mat size and get frobenus norm squared
    n_nodes = np.shape(M)[0]
    # f_nrm_sq = np.linalg.norm(M, ord = 'fro')**2      #unclear where used

    #ininitalize array of intervals[LEFT, RIGHT]
        #  NOTE: I think i'm doing
        #  this weird combo of matlab and python indexing. I'm going to continue 
        #  as is for now, using VECTOR INDEXING (x1, ..., x_nNodes)
    # Level j has exactly j blocks, 1 <= j <= n_nodes. Level 1 has only one block,
    # vector indices [1, n_nodes] (corresponding to interval [0,1]). 
    left = np.zeros((n_nodes, n_nodes))
    left[0,0] = 1
    right = np.zeros((n_nodes, n_nodes))
    right[0,0] = n_nodes
    # At each level n there are (n_nodes - n) potential "splits". First, the 
    # the algorithm finds all the potential splits at a given level, and follows
    # this by finding the soules basis vector associated with each split,
    # and then choose the soules basis vector which maximizes the value (see paper)
    # associated with a cluster containing the largest difference in edge connection

    #in moving from level n to n+1, we split exactly one block (re: support) into two

    #initialize output matrix
    Qn = np.zeros((n_nodes, n_nodes))
    Qn[0] = r1
    for n in range(n_nodes -1):     #gonna be a triple for loop ...
        #initialize arrays
        r = np.zeros((n_nodes-1, n_nodes)) #Set of possible wavelets at level n
        iblock = np.zeros(n_nodes - 1) # array of potential splitting indices
        coeff = np.zeros(n_nodes - 1) # coefficients for eeach possible wavelet
        where = np.zeros(n_nodes - 1) # location of splits

        #index of block overwhich creating a wavelet
        ib = 0      # indexed block? Potential error source
        for block in range(n+1):
            i0 = int(left[n,block])
            i1 = int(right[n,block])
            # we attempt to split the block at interval k such that
            #       [i0,i1] = [i0,k] U [k,i1]
            # if i0 == i1, that implies this is a leaf node
            if i0 < i1:
                for k in range(i0, i1):
                    #store the block index TODO: error check
                    iblock[ib] = block
                    #build the soules basis vector associated with this split
                    psi = build_vector(r1, i0, i1, k)
                    #pull the inner product coefficient
                    alpha = np.dot(psi, (M@psi))
                    #print(alpha)
                    #square said coefficient (for optimization)
                    coeff[ib] = alpha**2
                    #store the vector
                    r[ib] = psi
                    #record split location of vector
                    where[ib] = k
                    #incremement the wavelet index
                    ib += 1
        # compare the coefficients of the vectors across all blocks at level n
        # NOTE: this code sorts the ENTIRE array, and then pulls at random one of
        # the maximum values. Seems inefficient. Imma rewrite it instead :D
        # coeff_nonzero = coeff[np.nonzero(coeff)]
        # sorted_ind = np.argsort
        # Rewritten: Function to pull indices of only the maxes (within tolerance of 1E-9)
        # O(n), so faster than O(nlogn) of numpy sort?
        max_ind = util_funcs.all_max_ind(coeff, eps = 1E-9) 
        iw = np.random.choice(max_ind)
        # Get eigenvector at that location
        Qn[n+1] = r[iw]
        #update L and R ends of the blocks for level n+1 
            #  level n:
            #    [left(j)          k            right(j)] [left (j+1): right(j+1)]
            #    [left(j)        k][k+1         right(j)] [left (j+1): right(j+1)]
            # 
            #  level n+1:
            #    [left(j):right(j)][left(j+1):right(j+1)] [left (j+2): right(j+2)]

        # find index of block where split occured
        jk = int(iblock[iw])
        print(f"jk = {jk}\n")
        print(f"k = {where[iw]}\n")
        # for all the blocks before jk,the L and R indices remain the same
        if jk > 0:
            right[n+1, :jk] = right[n, :jk]
        left[n+1, :jk+1] = left[n, :jk+1]
        # split the jk-th block
        right[n+1, jk] = where[iw]
        left[n+1, jk+1] = where[iw] + 1
        right[n+1, jk+1] = right[n, jk]
        # shift the blocks after jk
        if jk < n:
            left[n+1, jk+2:n+2] = left[n, jk+1:n+1]
            right[n+1, jk+1:n+2] = right[n, jk:n+1]
    print(f"left = {left}\n")
    print(f"right = {right}\n\n")
    return Qn

A = genSBM(7,2,1,1)
# print(A)
test_output = topdown(A[0], 1/np.sqrt(7) * np.ones(7))
Qn = np.transpose(test_output)
test_id1 = Qn @ test_output
test_id2 = test_output @ Qn
print(Qn @ test_output)
print("debug")

def barycenter(AG, M):
    """
    _Description :    The function computes the Frechet mean of a sample of
                    T graphs, each of which of size n. The distance is the l2 norm
                    of the vector of the first M eigenvalues of the normalized graph
                    Laplacian.
    
    INPUT:
    
        AG: Tensor of n x n x T graphs of size n
    
        M: number of eigenvalues that are used to reconstruct the graph using the truncated Soules basis.
    
        n: size of the graphs
    
        T: number of graph in the sample
    
    OUTPUT:
    
        Qn: Soules basis; this is an orthonormal matrix of size n x n
    
        Bary: adjacency matrix of the barycentre graph reconstructed using the Soules basis
    
        Epd: indicator function of the block locations/nonzero entries in
            the barycentre graph; the diagonal has been removed.   
    
        EA: sample mean adjacency matrix; this is for comparison purposes
    
    _References :
    
    _Remarks : None
    
    _Author :                 Francois G. Meyer
    
    _Revisions History: 2025 Initial keying
    """
    T = np.shape(AG)

    return