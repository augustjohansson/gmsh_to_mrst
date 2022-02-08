function [G, h] = plot_gmsh(input, drawtags)

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    
    if nargin == 0
        input = 'Gdata.mat';
        drawtags = true;
    end

    if isa(input, 'struct')
        G = input;
    elseif isa(input, 'char')
        G = load(input);
    else
        error('Unknown input data of class %s', class(input));
    end

    G = computeGeometry(G, 'findNeighbors', true);

    figure, hold on
    if isfield(G.cells, 'tags')
        ut = unique(G.cells.tags);

        % Make sure background is yellow
        volumes = zeros(numel(ut), 1);
        for k = 1:numel(ut)
            idx = G.cells.tags == ut(k);
            volumes(k) = sum(G.cells.volumes(idx));
        end
        [volumes, idx] = sort(volumes, 'descend');
        max_vol_tag = idx(1);
        
        colors = [lines(7); rand(numel(ut)-7, 3)];
        
        for k = 1:numel(ut)
            if ut(k) == max_vol_tag
                color = [1,1,0]; % yellow
            else
                color = colors(k,:);
            end
            idx = G.cells.tags == ut(k);
            dispif(mrstVerbose, 'Cell tag %d: %d cells\n', ut(k), sum(idx))
            h{1} = plotGrid(G, idx, 'facecolor', color);
        end
    else
        h{1} = plotGrid(G);
    end

    if drawtags
        if isfield(G.faces, 'tags')
            ut = unique(G.faces.tags);
            default_face_tag = 0;
            dispif(mrstVerbose, 'Ignoring face tag %d\n', default_face_tag);
            ut(ut==default_face_tag) = [];
            colors = [lines(7); rand(numel(ut)-7, 3)];
            for k = 1:numel(ut)
                idx = G.faces.tags == ut(k);
                dispif(mrstVerbose, 'Face tag %d: %d faces\n', ut(k), sum(idx))
                h{2} = plotFaces(G, idx, 'linewidth', 2, 'edgecolor', colors(k,:));
            end
        end
    end

end


