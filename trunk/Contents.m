% Toolbox MatrixClust
% Copyright (C) 2008 Tom Michoel
% Version 2008-08-04
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
%% This set of programs is free software; you can redistribute it and/or 
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation; either version 2 of the
%% License, or (at your option) any later version.
%%
%% This set of programs is distributed in the hope that it will be
%% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License along
%% with this program; if not, write to the Free Software Foundation, Inc.,
%% 675 Mass Ave, Cambridge, MA 02139, USA.
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
%% You are welcome to use the code for your research under the terms of the
%% license. However, please acknowledge its use with the following
%% citation:
%%
%% A. Joshi, Y. Van de Peer and T. Michoel; Analysis of a Gibbs sampler
%%   method for model based clustering of gene expression data;
%%   Bioinformatics 24(2):176-183 (2008).
%%
%% Author contact information:
%%
%%  Tom Michoel
%%
%%  Bioinformatics & Evolutionary Genomics
%%  Department of Plant Systems Biology
%%  VIB, Ghent University
%%  Technologiepark 927
%%  B-9000 Gent, Belgium
%%
%%  Email: tom.michoel@psb.ugent.be
%%  URL:   http://www.psb.ugent.be/~tomic/
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% MatrixClust is a Matlab toolbox for fuzzy clustering of a symmetric
% matrix, typically the weighted adjacency matrix of an undirected network.
%
% 1. INSTALLATION
% ================
%
% To install, unpack in a directory DIR and add to Matlab search path or
% change into DIR and run: 
%
% >> Install
%
%
% 2. FUNCTIONS
% =============
%
% matrixClustSym :
%   
%   Improved version of the fuzzy clustering algorithm described in
%   [A. Joshi, Y. Van de Peer and T. Michoel; Bioinformatics 24:176 (2008)]
%
%   For usage information, type:
%
%   >> doc matrixClustSym
%
% pcutoff :
%   
%   Sets all elements of a matrix below certain cutoff to zero.
%
%   For usage information, type:
%
%   >> doc pcutoff
%
% optCutoff :
%
%   Determines optimal probability cutoff in fuzzy clusters by maximizing
%   the ratio of the (weighted) number of edges to the (weighted) number
%   of nodes.
%
%   For usage information, type:
%
%   >> doc optCutoff
%