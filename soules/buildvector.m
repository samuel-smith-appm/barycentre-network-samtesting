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
%%_Module_Name : buildvector.m
%%
%%_Description : construct a Soules vector. 
%%              
%%
%%_References :
%%
%%_Remarks : None
%%
%%_Author :                 Francois G. Meyer
%%
%%_Revisions History: 2025 Initial keying
%%
%%==============================================================================


function [v] = buildvector (x,i0,i1,k)

  %% input: x vector of length sizex

  %% [i0,i1] set of indices of the block from which we construct the basis vector
  %% k: index where we split: [i0,k] on one side and [k+1,i1] on the other side

  %% output: the basis vector corresponding to the node

  if (k >= i1)
    fprintf('\n\t function buildvector: runtime error');
    fprintf('\n\t index of the split is too large %d; last block index is %d', k, i1);
    return;
  end
  
  node  = x(i0:i1,1);

  %%        +--------+  
  %% left:  |        |
  %%        |        +---------------------
  %%        i0       k                   i1
  %%                 +--------------------+
  %% right:          |                    |
  %%        ---------+                    |

  left  = cat (1, x(i0:k,1), zeros (i1 -k,1));
  right = cat (1, zeros (k-i0 + 1, 1), x(k+1:i1,1));
  
  n = norm (node);

  nl = norm (left);
  nr = norm (right);

  v = zeros (size(x));

  %%        +--------+  
  %%        |        |
  %%  v: __ i0 _____ k __________________ i1 ___
  %%                 |                    |
  %%                 +--------------------+

  v(i0:i1) = (nr/nl.*left - nl/nr.*right)./n;

  return;
end

