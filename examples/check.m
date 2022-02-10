clear all
close all

models = {'tets', 'hexes', 'prisms'};

for k = 1:numel(models)
    models{k}
    pyname = sprintf('%s.py', models{k});
    % pyrunfile(pyname);
    [status, cmdout] = system(sprintf('python %s', pyname));
    assert(status == 0, [num2str(status), cmdout]);

    mname = sprintf('%s.m', models{k});
    G = gmsh_to_mrst(mname);
    assert(all(G.cells.volumes > 0));
    assert(all(G.faces.areas > 0));

end
