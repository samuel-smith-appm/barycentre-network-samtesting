%%==============================================================================
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
%%_Module_Name : IHfast
%%
%%_Description : construct an inhomogeneous random graph based on independent
%% Bernoulli symmetric matrix with chances p_ij  
%%
%%
%%_References :
%%
%%_Remarks : This is fast.
%%
%%_Author :                 Francois G. Meyer
%%
%%_Revisions History: 2025 Initial keying
%%
%%==============================================================================
%%

function [A] = IHfast (P)

  petit = eps(1);
  
  nVertices = size(P,1);
  
  A = zeros (nVertices, nVertices);
  
  %% take upper diagonal without the diagonal
  
  topindices =  triu(true(size(A)), 1)
  
  %% generate independent Bernoulli entries until we have a connected graph
  %% we do not test connectivity, we simply check no isolated points.
  %% In theory this is the same with high probability.
  P(topindices)
  code = -1;
  while (code == -1)
    
    topA = lebernoulliIH (P(topindices))
    
    A (topindices) = topA;
    A = A';			%% transpose the matrix; 
    
    A (topindices) = topA;	%% mirror symmetric copy
    
    dA = A*ones(nVertices,1);

    if (min (abs(dA)) < petit)
      fprintf('\nWarning the graph has isolated nodes');
    else
      code = 1;
    end
end

