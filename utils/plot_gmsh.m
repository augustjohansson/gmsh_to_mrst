function [G, h] = plot_gmsh(filename)

    if nargin == 0
        filename = 'Gdata.mat';
    end

    G = load(filename);
    G = computeGeometry(G, 'findNeighbors', true);

    figure, hold on
    if isfield(G.cells, 'tags')
        ut = unique(G.cells.tags);
        colors = [lines(7); rand(numel(ut)-7, 3)];
        
        for k = 1:numel(ut)
            idx = G.cells.tags == ut(k);
            dispif(mrstVerbose, 'Cell tag %d: %d cells\n', ut(k), sum(idx))
            h{1} = plotGrid(G, idx, 'facecolor', colors(k,:));
        end
    else
        h{1} = plotGrid(G);
    end

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


