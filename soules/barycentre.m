%%=================================================================================================
%%
%% Copyright (C) 2025 Francois G. Meyer <FMeyer@Colorado.Edu>
%%
%% All rights reserved.
%%
%% Redistribution and use in source and binary forms, with or without
%% modification, are permitted provided that the following conditions are met:
%% 
%%   a. Redistributions of source code must retain the above copyright notice,
%%      this list of conditions and the following disclaimer.
%% 
%%   b. Redistributions in binary form must reproduce the above copyright
%%      notice, this list of conditions and the following disclaimer in the
%%      documentation and/or other materials provided with the distribution.
%% 
%%   c. Neither the name of the copyright holders nor the names of any
%%      contributors to this software may be used to endorse or promote products
%%      derived from this software without specific prior written permission.
%% 
%% 
%% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS on an
%% "AS IS" basis. THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO
%% REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  BY WAY OF EXAMPLE, BUT
%% NOT LIMITATION, THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO AND
%% DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR
%% ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE
%% ANY THIRD PARTY RIGHTS.
%% 
%% THE COPYRIGHT HOLDERS AND CONTRIBUTORS SHALL NOT BE LIABLE TO LICENSEE OR
%% ANY OTHER USERS OF THIS SOFTWARE FOR ANY INCIDENTAL, SPECIAL, OR
%% CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR
%% OTHERWISE USING THIS SOFTWARE, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED
%% OF THE POSSIBILITY THEREOF.
%%
%%
%%_Module_Name : barycentre.m
%%
%%_Description :    The function computes the Frechet mean of a sample of
%%                  T graphs, each of which of size n. The distance is the l2 norm
%%                  of the vector of the first M eigenvalues of the normalized graph
%%                  Laplacian.
%%
%%  INPUT:
%%
%%    AG: Tensor of n x n x T graphs of size n
%%
%%    M: number of eigenvalues that are used to reconstruct the graph using the truncated Soules basis.
%%
%%    n: size of the graphs
%%
%%    T: number of graph in the sample
%%
%%  OUTPUT:
%%
%%    Qn: Soules basis; this is an orthonormal matrix of size n x n
%%
%%    Bary: adjacency matrix of the barycentre graph reconstructed using the Soules basis
%%
%%    Epd: indicator function of the block locations/nonzero entries in
%%         the barycentre graph; the diagonal has been removed.   
%%
%%    EA: sample mean adjacency matrix; this is for comparison purposes
%%
%%_References :
%%
%%_Remarks : None
%%
%%_Author :                 Francois G. Meyer
%%
%%_Revisions History: 2025 Initial keying
%%
%%=================================================================================================

function [Qn, Bary, Epd, EA] = barycentre (AG,M,n,T)
  %%
  %%
  %%
  
  DISPLAY = 0;

  %%
  %% Sample mean adjacency matrix; we need it for the Soules basis
  %%

  EA = mean (AG,3);

  %%
  %% display the chance matrix for the SBM, and the SBM
  %%

  if (DISPLAY)
    f=figure;  f.Theme = "light";
    imagesc(EA); axis square, axis off; colormap ('jet'); colorbar;
    fig2pdf (f, 'sample mean A');
  end

  oneshot = (T == 1);

  if (oneshot)
    K = AG;

    dK   = sum(K,2);		%% compute the degree vector
    iDK2 = diag(1./sqrt(dK));	%% normalization matrix = 1./sqrt(degree matrix)
    
    Kn = iDK2 * K * iDK2;	
    Kn = 0.5.*(Kn + Kn');	%% make it symmetric, for good measure
    
    %% Compute the eigenvalues of \hat{K}

    lambda = 1 - eig(Kn);								

    %% sort the eigennvalues by increasing order, and get rid of the bulk eigenvalues

    eigL  = sort(lambda,'ascend');
    eigL2 = eigL (1:M);

  else
    %% compute the  eigenvalues of each normalized adjacency matrix
    
    lamb = zeros (n,T);
    Kn   = zeros (n,n);
    K    = zeros (n,n);
    
    for t=1:T
      K = AG(:,:,t);
      
      dK   = sum(K,2);
      iDK2 = diag(1./sqrt(dK));
      
      Kn = iDK2 * K * iDK2;
      Kn = 0.5.*(Kn + Kn');
      
      lamb (:,t) = eig(Kn);
      
    end
    
    %% this is the mean spectrum
    
    lambda = 1 - mean (lamb,2);
    
    %% sort the eigenvalues by increasing order, and get rid of the bulk eigenvalues
    
    eigL  = sort(lambda,'ascend');
    eigL2 = eigL (1:M);

    if (DISPLAY)
      nbins = 200;
      f= figure;histogram (eigL, nbins, 'Normalization','probability');
      fig2pdf (f,'distribution_eig_L');
    end
    
  end

  %%==================================================================================================
  %%
  %% Compute the best Soules basis
  %%
  %%==================================================================================================

  %% first eigenvector of the normalized graph Laplacian

  x = (1/sqrt(n)).*ones (n,1);

  %%
  %% allocate memory
  %%
  Qn = zeros (n, n);

  %%
  %% build the Soules basis by exploring the vector of cuts
  %%

  Qn = topdown (EA,x);

  %%==================================================================================================
  %%
  %% remove the columns of the Soules basis associated to these eigenvalues  and compute EM (see paper)
  %% We also threshold EM to get an estimate of the blocks geometry, EP.
  %% Finally, we remove the diagonal in EP, and rerurn Epd.
  %%==================================================================================================

  Q = Qn (:,1:M);
  EM  = Q*Q';
  EP  = EM > 1e-6;
  Epd = EP - diag (diag(EP));

  clear ('EP');
  
  if (DISPLAY)
    figure;imagesc(EM);title ('EM'); colormap('jet'); colorbar;
  end

  %%=================================================================================================
  %%
  %% estimate of the degree
  %% the key idea is to compute the average of the degrees in each block
  %%
  %%=================================================================================================

  d1     = sum (EA,2);
  blcksz = sum(Epd,2);
  dbar   = Epd*d1;
  dbar   = dbar./blcksz;

  Dbar   = diag(dbar);
  D2bar  = diag(sqrt(dbar));

  %%=================================================================================================
  %%
  %% test the reconstruction of the normalized Laplacian using the
  %% first M Soules eigenvectors where M is the number of communities
  %% reconstruct an adjacency matrix using the expected degree matrix 
  %%
  %%=================================================================================================

  L2 =  (Q * diag(eigL2) *Q');
  A2 =  D2bar * (EM - L2) * D2bar;

  A2 = 0.5.*(A2 + A2');
  A2 = A2 - diag(diag(A2));

  Atrunc = A2;

  clear ('A2', 'Dbar',  'EM');

  %%=================================================================================================
  %%
  %% reconstruction using all the eigenvectors
  %% The following code works with precision 1e-15 when A is replaced with its expectation: P.
  %% To wit sup|A-AT| and sup |A-AL| are <= 1e-15
  %%
  %%=================================================================================================

  In = eye(n, n);

  L3 = Qn * diag(eigL) *Qn';
  A3 = (D2bar * (In - L3) * D2bar);

  A3 = 0.5.*(A3 + A3');
  A3 = A3 - diag(diag(A3));

  Afull = A3;

  clear ('A3', 'D2bar', 'L3');

  %% if (DISPLAY)
  %%   f = figure;imagesc (A3); colorbar; colormap('jet');
  %%   axis square, axis off;
  %%   fig2pdf (f,'full reconstruction');
  %% end

  %%
  %% we return the truncated reconstruction
  %%

  Bary = Atrunc;

  %%
  %% we return the full reconstruction
  %%
  
  %% Bary = 0.5*(Afull + Afull');

end





