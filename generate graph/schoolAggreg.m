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
%_Module_Name: schoolAggreg.m
%
%_Description: load the unweighted aggregated two-day school dataset
%
%   Input parameters: 
%          min0: in minutes, time index of the first network to analyze.
%          min1: in minutes, time index of the last network to analyze.
%          aggMin: length of the time window over which we aggregate the graphs,
%          measured in minutes. See [1,2,] for in-depth analysis of the optimal choice
%          for this parameter. 
%
%   Output parameters:
%
%         G : tensor of size n x n x T of T adjacency matrices.
%         leseleves: (for display purposes) class label of each student
%         miniP: average edge probability for each community, and accross 
%         communities (last entry).  
%
%_References:
%
%  [1] Bovet, A., Delvenne, J.C. and Lambiotte, R., 2022. Flow
%      stability for dynamic community detection. Science advances, 8(19),
%      p.eabj3063
%
%  [2] Failla, A., Cazabet, R., Rossetti, G. and Citraro, S.,
%      2024. Describing group evolution in temporal data using
%      multi-faceted events. Machine Learning, 113(10), pp.7591-7615. 
%
%_Remarks: None
%
%_Author:                 Francois G. Meyer
%
%_Revisions History: 2025 Initial keying
%
%==============================================================================
function [G, leseleves, miniP] = schoolAggreg (min0,min1,aggMin)

  DISPLAY = 0;

  name = sprintf ('every_%03d_%03d_%03d.mat', aggMin, min0, min1);
  fileName = sprintf('./data/%s', name);

  load (fileName, 'leclock', 'U', 'B', 'leseleves', 'clssSz');

  M = 10;			% number of communities

  G = U;			% tensor of dynamic networks
  
  A = mean (G,3);		% mean (over time) network
  
  n = size (A,1);		% number of students

  One  = ones (n,n);

  if (DISPLAY)
    figure;
    drawgraph (A, leseleves);
  end

  %=================================================================================================
  %
  % statistics on the primary school dataset;
  % we compute:
  %             - the edge probability within each class
  %             - the (overal) edge probability between classes
  %
  %=================================================================================================

  aved = zeros (M,1);
  
  outEdges = 0;
  outAll   = 0;
  
  deb = 1;
  
  for (m=1:M)
    b = clssSz (m);		                  % size of current block
    fin = deb + b - 1;		                  % indices of the block vertices in P
    idx = [deb:fin]; 

    outAll = outAll + (n - b) * b;                % number of all possible edges from the block to the rest of the world

    inEdges  = (sum(A(idx, idx),"all"));          % total number of edges inside the block
    aved (m) = inEdges / ((b-1)*(b-1));	          % average degree in the block

    allEdges = (sum(A(idx,:),"all"));             % all edges from the block
    outEdges = outEdges + allEdges - inEdges;     % actual edges from the block to the rest of the world
    
    deb = fin + 1;		                  % next block
  end

  % miniP (1:M) = average edge probability inside each community
  % miniP (M+1) = average edge probability across communities
  
  miniP = zeros (M+1,1);
  miniP (1:M) = aved;
  
  q0 = outEdges/outAll;
  miniP (M+1) = q0;

end


