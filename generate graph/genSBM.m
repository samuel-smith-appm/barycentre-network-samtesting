%______________________________________________________________________________
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
%_Module_Name :
%
%_Description : Generate one random sample from an SBM model.
%               The function handles two geometries:
%               1) an unbalanced model with 4 blocks: this model is trigerred by equalSize = 0
%               2) a balancd model of M blocks; this model is trigerred by equalSize = 1
%               All models can have arbitrary sizes (number of nodes).
%               The edge probability is of the form O(log^2 (n)/n) for the edge probability in a 
%               community, and O(log(n)/n) for edges between communities.
%
%               Input:
%                      n: graph size
%                      M: number of blocks
%                      equalSize: 1 if balanced SBM
%                      equalChance: 1 if same edge probability in all the blocks
%
%               Output:
%                      P: edge probability matrix: n x n
%                      p: =[p1,..,pM,q]
%                      sizeofBlock = n/M
%
%_References :
%
%_Remarks : None
%
%_Author :                 Francois G. Meyer
%
%_Revisions History: 2025 Initial keying
%
%______________________________________________________________________________



function  [P, p, sizeofBlock] = genSBM (n, equalSize, equalChance, M)

  %
  % M is the number of blocks; this is only relevant if equalSize = 1 and equalChance = 1.
  %
  
  addpath ('./..');

  switch (equalSize)
    case 0

      %number of communities

      if (M ~= 4)
	fprintf ('\n Warning: you requested %d blocks, but the unbalanced model has only 4 blocks\n', M);
	M = 4;
      end

      %  communities have different sizes

      chunk = round(n/24);
      sizeofBlock (1) = 3*chunk;
      sizeofBlock (2) = 7*chunk;
      sizeofBlock (3) = 5*chunk;
      sizeofBlock (4) = n - sum(sizeofBlock(1:3));

    case 1
      % This choice of communities is designed to test the estimates provided in lowe 2024
      % These estimates are rather crude; they provide one digit precision

      % all blocks have the same size

      blockSize     = floor (n/M);
      sizeLastBlock = n - (M-1)*blockSize;
      sizeofBlock   = blockSize * ones(M,1);
      sizeofBlock (end) = sizeLastBlock;
      
  end

  % the coefficients below guarantee a good separation between the M dominant eigenvalues of the normalized 
  % adjacency matrix and the bulk; see lowe 2024 for details.

  if (equalChance)
    chance = 3;
    coeff  = chance *ones (M,1);
  else
    c     = linspace (3,M+2,M)';	% 1,2,..,M
    coeff = c(randperm (M));	        % randomly permute these coefficients
  end
  cq = 1;

  % load the edge probability within (p) and across (q) communties
  
  p      = ((log(n).^2)/n) .* ones (M+1,1);
  p      = coeff.* p (1:M);
  p(M+1) = cq  * log(n)/n;

  % fill in the chance matrix P
  % one block at aa time

  P =  p (M+1) * ones (n,n);

  deb = 1;
  for (m=1:M)
    fin = deb + sizeofBlock(m) - 1;
    idx = [deb:fin]; 
    P (idx, idx) = p(m);
    deb = fin + 1;
  end

  % there are no self loop, so we remove the diagonal
  
  P = P - diag(diag(P));

end
