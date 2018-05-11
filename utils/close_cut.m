function [Torus_ori] = close_cut(Torus)
%
cId = 0;
numP = size(Torus.vertexPoss, 2);
IDX = zeros(1, numP);
for pId = 1 : numP
    if IDX(pId) == 0
        cId = cId + 1;
        IDX(pId) = cId;
        for j = (pId+1) : numP
            if IDX(j) == 1
                continue;
            end
            d = Torus.vertexPoss(:, pId) - Torus.vertexPoss(:, j);
            if d'*d < 1e-10
                IDX(j) = cId;
            end
        end
    else
        continue;
    end
end
Torus_ori.faceVIds = IDX(Torus.faceVIds);
Torus_ori.vertexPoss = zeros(3, cId);
for pId = 1 : numP
    Torus_ori.vertexPoss(:, IDX(pId)) = Torus.vertexPoss(:, pId);
end
for i = 1:cId
    count(i) = length(find(IDX ==i));
end
%
clusIds = find(count == 2);
G = mesh_graph(Torus);
%
for i = 1 : length(clusIds)
    vids = find(IDX == clusIds(i));
    sId = vids(1);
    tId = vids(2);
    [i1, i2, i3] = graphshortestpath(G, sId, tId);
    Pairs(1,i) = sId;
    Pairs(2,i) = tId;
    Pairs(3,i) = i1;
end
h = 10;