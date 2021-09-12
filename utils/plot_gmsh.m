function plot_gmsh(filename, titletag)

    if nargin == 0
        filename = 'Gdata.mat';
        titletag = '';
    end
    if nargin == 1
        titletag = '';
    end
    
    mrstVerbose on

    G = load(filename);
    G = computeGeometry(G, 'findNeighbors', true);

    figure, hold on
    if isfield(G.cells, 'tags')
        ut = unique(G.cells.tags);
        colors = [lines(7); rand(numel(ut)-7, 3)];
        
        for k = 1:numel(ut)
            idx = G.cells.tags == ut(k);
            dispif(mrstVerbose, 'Cell tag %d: %d cells\n', ut(k), sum(idx))
            plotGrid(G, idx, 'facecolor', colors(k,:));
        end
    else
        plotGrid(G)
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
            plotFaces(G, idx, 'linewidth', 2, 'edgecolor', colors(k,:));
        end
    end

    grid on
    axis equal tight
    %title(sprintf('%s cells %d, faces %d, nodes %d', titletag, G.cells.num, G.faces.num, G.nodes.num))

end

