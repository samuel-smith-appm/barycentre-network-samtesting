%==============================================================================
%
% Copyright (C) 2025 Francois G. Meyer <FMeyer@Colorado.Edu>
%
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
%   a. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
% 
%   b. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
% 
%   c. Neither the name of the copyright holders nor the names of any
%      contributors to this software may be used to endorse or promote products
%      derived from this software without specific prior written permission.
% 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS on an
% "AS IS" basis. THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO
% REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  BY WAY OF EXAMPLE, BUT
% NOT LIMITATION, THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO AND
% DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR
% ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE
% ANY THIRD PARTY RIGHTS.
% 
% THE COPYRIGHT HOLDERS AND CONTRIBUTORS SHALL NOT BE LIABLE TO LICENSEE OR
% ANY OTHER USERS OF THIS SOFTWARE FOR ANY INCIDENTAL, SPECIAL, OR
% CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR
% OTHERWISE USING THIS SOFTWARE, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED
% OF THE POSSIBILITY THEREOF.
%
%
%_Module_Name : topdown
%
%_Description : explores the binary tree of Soules bases, from the top (coarsest)
%               level, down to the finest level. The function computes the best Soules
%               basis , which provides the sparest expansion of the matrix A.
%
%               The core function is buildvector that constructs a new Soules vector, based
%               on entries in r1, and the geometry in the Soules tree that is being split.
%
%               Input:
%                       A the sample mean adjacency matrix, of size nNodes x nNodes
%                       r1 first Soules vector: nonnegative vector; typically r1 = n^{-1/2} one
%               Output:
%                       Qn  orthogonal matrix to store all the basis vectors 
%                           only nNodes - 1 columns are non zero in Qn; the 
%
%_References :
%
%_Remarks : None
%
%_Author :                 Francois G. Meyer
%
%_Revisions History: 2025 Initial keying
%
%==============================================================================

function  [Qn] = topdown  (M,r1)
  
  % size of the adjacency matrix
  nNodes = size(M,1);								

  % Frobenius norm
  norm2 = norm (M,"fro").^2;							
  
  left  = zeros (nNodes,nNodes);
  right = zeros (nNodes,nNodes);
  
  % initialize the array of intervals [left,right] 
  %
  %   Level j  has exactly j blocks, 1 <= j <= nNodes
  %
  %   Each block j has the following boundary elements:
  %
  %   [left(j):right(j)][left(j+1):right(j+1)]
  %      block j           block j + 1 

  % at the top level, there is only one block, [1,nNodes]
  left  (1,1) = 1;
  right (1,1) = nNodes;		                                                

  %
  % we explore top down: going down the levels
  % to go from level n to level n+1 we need split exactly one block
  % we choose the block which yields the largest wavelet coefficient
  %
  % This version of the code allocates the memory beforehand

  % matrix of Soules vectors
  Qn  = zeros (nNodes, nNodes);
  Qn(:,1) = r1;

  % in principle there are nNodes coefficients, but we need to remove one 
  % coefficient per block, so we really need nNodes -n coefficients at level n

  % set of all possible wavelets at level n
  r      = zeros (nNodes, nNodes - 1);

  % for a given level n, index of the block wherein we split
  iblock = zeros (nNodes - 1, 1);

  % coefficients for every possible wavelet
  coef   = zeros (nNodes - 1, 1);

  % absolute location of where the split occurs in the block at level n
  where  = zeros (nNodes - 1, 1);

  for n = 1:nNodes-1                                
    % reset the arrays

    r(:)   = 0;

    iblock (:) = 0;
    coef (:)   = 0;
    where (:)  = 0;
    
    % index of the block wherein we construct a wavelet
    ib = 1;

    %
    % For each block bl at level n, we try to split the block using an index k 
    % in [left(bl),..,right(bl)-1] = [i0,i1-1]
    % We construct the eigenvector associated with the split, and compute the
    % correlation with the adjacency matrix
    %
    
    % for each block
    for bl = 1:n								

      % begining and end of the block
      
      i0 = left  (n, bl);							
      i1 = right (n, bl);							

      % we try to split the block at  k, i0 <= k < i1
      % [i0,i1] = [i0,k] U [k+1,i1] 
      % if i0 == i1 then the block is a leaf and we cannot create a new eigenvector
      
      if (i0 < i1)
	% this block is not a leaf, so we try to split, wherever we can

        for k = i0:i1-1	

          % for the block index ib, store the index of the block
          iblock(ib) = bl; 
          
          % this vector oscillates over the block: psi > 0 on [i0,k], and psi < 0 on [k+1,i1]
          psi = buildvector (r1, i0, i1, k); 

                % coefficient associated with the tensor product of the Soules vector

          alpha = psi' * M * psi;
          coef(ib)  = alpha*alpha;

          % store the vector
          r (:,ib) = psi;

          % absolute location of the split in the block bl
          where (ib) = k;

          % increment the wavelet index
          ib = ib+1;
        end
      end
    end
    
    % compare all vectors, across all blocks at level n
    [tmp,idx] = sort (coef (1:ib-1),'descend');

    % the following three lines handle the case where there are multiple 
    % maximum coefficients that are equal (whithin 10-9)
    % pick one vector randomly, if they are all equal

    lesindex = idx (find (abs(tmp - tmp(1)) < 1e-9));
    if (length(lesindex) == 1)
      iw = lesindex;
    else 
      % pick one of the maxima at random
      iw = randsample (lesindex,1);
    end

    % fprintf ('\t level = %d\t cut = %d\n', n, where(iw));

    % get the eigenvector at that location 
    v = r (:,iw);
    Qn(:,n+1) = v;

    % update the left and right ends of the blocks for the next level: n+1
    %
    % level n:
    %   [left(j)          k            right(j)] [left (j+1): right(j+1)]
    %   [left(j)        k][k+1         right(j)] [left (j+1): right(j+1)]
    %
    % level n+1:
    %   [left(j):right(j)][left(j+1):right(j+1)] [left (j+2): right(j+2)]
    
    % find the index of the block where the split happened
    
    jk = iblock (iw);

    %
    % for all the blocks before jk, the left and right indices remain the same.
    %
    if (jk > 1)
      right (n+1, 1:jk-1) = right (n,1:jk-1);
    end
    left  (n+1, 1:jk)   = left  (n,1:jk);

    % split the block jk

    right (n+1, jk)     = where (iw);
    left  (n+1, jk + 1) = where (iw) + 1;
    right (n+1, jk + 1) = right (n, jk);

    %
    % for all the blocks after jk, the left and right indices remain the same.
    %
    if (jk < n)
      left  (n+1, jk+2:n+2) = left  (n, jk+1:n+1);
      right (n+1, jk+2:n+2) = right (n, jk+1:n+1);
    end
    
  end				% next level down


  return;

