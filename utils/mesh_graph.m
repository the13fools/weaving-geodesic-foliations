function [G] = mesh_graph(Torus)
numP = size(Torus.vertexPoss, 2);
G = sparse(numP, numP);
for fId = 1 : size(Torus.faceVIds, 2)
    v1Id = Torus.faceVIds(1, fId);
    v2Id = Torus.faceVIds(2, fId);
    v3Id = Torus.faceVIds(3, fId);
    G(v1Id, v2Id) = 1;
    G(v2Id, v1Id) = 1;
    G(v1Id, v3Id) = 1;
    G(v3Id, v1Id) = 1;
    G(v2Id, v3Id) = 1;
    G(v3Id, v2Id) = 1;
end