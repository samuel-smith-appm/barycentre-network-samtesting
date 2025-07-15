"""
______________________________________________________________________________

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
import matplotlib.pyplot as plt


def genSBM(n, M, equal_size = True, equal_chance = True):
    """
    _Module_Name :

    _Description : Generate one random sample from an SBM model.
                The function handles two geometries:
                1) an unbalanced model with 4 blocks: this model is trigerred by equalSize = 0
                2) a balancd model of M blocks; this model is trigerred by equalSize = 1
                All models can have arbitrary sizes (number of nodes).
                The edge probability is of the form O(log^2 (n)/n) for the edge probability in a 
                community, and O(log(n)/n) for edges between communities.

                Input:
                        n: graph size
                        M: number of blocks
                        equalSize: 1 if balanced SBM
                        equalChance: 1 if same edge probability in all the blocks

                Output:
                        P: edge probability matrix: n x n
                        p: =[p1,..,pM,q]
                        sizeofBlock = n/M

    _References :

    _Remarks : None

    _Author :                 Francois G. Meyer

    _Revisions History: 2025 Initial keying
"""
    if equal_size:
        #This choice of communities is designed to test the estimates provided in lowe 2024
        #These estimates are rather crude; they provide one digit precision

        # observe that we do not enforce that n/M is an integer, and instead we simply
        # truncate the last block
        block_size = int(np.floor(n/M))
        final_block_size = n - (M-1)*block_size
        size_of_block = np.ones(M) * block_size
        size_of_block[-1] = final_block_size
    else:
        # this part is confusing? seems like equal_size = False only allows for this very
        # specific SBM model
        if M != 4:
            print(f"Warning: requested {M} blocks, but unbalanced model is only supported for 4 block!\n")
        M=4
        size_of_block = np.zeros(4)
        chunk = np.round(n/24)
        size_of_block[0] = 3*chunk
        size_of_block[1] = 7*chunk
        size_of_block[2] = 5*chunk
        size_of_block[3] = n + (chunk * 15)
    # the coefficients below guarantee a good separation between the M dominant eigenvalues of the normalized 
    # adjacency matrix and the bulk; see lowe 2024 for detail 
    if equal_chance:
        chance = 3
        coeff = chance * np.ones(M)
    else:
        coeff = np.linspace(3, M+2, M)
        coeff = np.random.permutation(coeff)
    cq = 1
    # load edge probability p within communities and q between communities
    p = ((np.log(n)**2)/n) * np.ones(M) * coeff
    p = np.append(p, np.log(n)/n)
    # initialize probability matrix P
    P = p[-1] * np.ones((n,n))
    # fill in block by block
    init_ran = 0
    for j in range(M):
        fin_ran = int(init_ran + size_of_block[j])
        P[init_ran:fin_ran, init_ran:fin_ran] = p[j]
        init_ran = fin_ran
    #remove diagonal
    P = P - np.diag(np.diag(P))



    return P, p, size_of_block

def lebernoulliTH(P):
    """
    _Module_Name : lebernoulliH
    
    _Description : generates an n x m matrix with Bernoulli entries 
    
    
    _References : 
    
    _Remarks : This is faster than the MATLAB version.
    
    _Author :                 Francois G. Meyer
    
    _Revisions History: 2025 Initial keying
    """
    [n,m] = np.shape(P)
    A = np.where(np.random.rand(n,m) < P, 1, 0)
    return A
def IHfast(P):
    """
    _Module_Name : IHfast
    
    _Description : construct an inhomogeneous random graph based on independent
     Bernoulli symmetric matrix with chances p_ij  
    
    
    _References :
    
    _Remarks : This is fast.
    
    _Author :                 Francois G. Meyer
    
    _Revisions History: 2025 Initial keying
    """
    n_nodes = np.shape(P)[0]
    A = np.zeros((n_nodes, n_nodes))
    top_indices = np.where(np.triu(np.ones((n_nodes, n_nodes)), k=1)==1, 1, 0)
    #  generate independent Bernoulli entries until we have a connected graph
    #  we do not test connectivity, we simply check no isolated points.
    #  In theory this is the same with high probability.
    code = -1
    while code == -1:
        topA = lebernoulliTH(P[top_indices])
        code = 1

    return

P = 0.5* np.ones((6,6))
IHfast(P)
# A = lebernoulliTH(P)
# print(P)
# print(A)


# test = False
# if test:

#     p1 = genSBM(25, 5, 1, 0)[0]
#     p2 = genSBM(25, 6, 1, 1)[0]
#     p3 = genSBM(25, 4, 0, 1)[0]
#     p4 = genSBM(25, 4, 0, 0)[0]

#     plt.imshow(p1)
#     plt.colorbar()
#     plt.show()
#     plt.imshow(p2)
#     plt.colorbar()
#     plt.show()
#     plt.imshow(p3)
#     plt.colorbar()
#     plt.show()
#     plt.imshow(p4)
#     plt.colorbar()
#     plt.show()
