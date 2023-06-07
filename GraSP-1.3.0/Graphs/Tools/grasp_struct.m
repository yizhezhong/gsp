%Defines and returns an empty graph structure.
%
%   GRASP_STRUCT()
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-2016)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018-2019)
% 
% benjamin.girault@ens-lyon.fr
% benjamin.girault@usc.edu
% 
% This software is a computer program whose purpose is to provide a Matlab
% / Octave toolbox for handling and displaying graph signals.
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,  more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.

function g = grasp_struct
    g.A = 0;                         % Adjacency matrix. g.A(j,i) is the weight of the edge from i to j
    g.A_layout = 0;                  % Adjacency matrix to plot
    g.layout = 0;                    % Graph layout (2D or 3D)
    g.distances = 0;                 % Matrix of distances between any pair of nodes
    g.node_names = {};               % Cell of strings for the names of nodes
    g.M = 0;                         % Matrix of the graph variation operator (e.g. the Laplacian)
    g.Q = 0;                         % Matrix of the graph signal inner product (e.g. identity)
    g.Z = 0;                         % Matrix of the fundamental matrix of the graph
    g.fourier_version = 'n.a.';      % Version of the graph Fourier transform used to get eigvals, F and Finv
    g.eigvals = 0;                   % Graph frequencies of the graph
    g.F = 0;                         % Fourier matrix
    g.Finv = 0;                      % Inverse Fourier matrix
    g.T = 0;                         % Translation operator
    g.background = '';               % Background image file path
    g.show_graph_options = struct(); % Default options for grasp_show_graph
end
