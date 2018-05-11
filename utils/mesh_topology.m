function [edges, adjFaces, faceEIds] = mesh_topology(Mesh)
% Compute the topology of the triangular mesh
%
numV = size(Mesh.vertexPoss, 2);
numF = size(Mesh.faceVIds, 2);
v1Ids = Mesh.faceVIds(1, :);
v2Ids = Mesh.faceVIds(2, :);
v3Ids = Mesh.faceVIds(3, :);
%
rows_A_vv = [v1Ids, v2Ids, v3Ids];
cols_A_vv = [v2Ids, v3Ids, v1Ids];
vals_A_vv = ones(1, 3*numF);
A_vv = sparse(rows_A_vv, cols_A_vv, vals_A_vv, numV, numV);
A_vv = A_vv + A_vv';
[rows, cols, vals] = find(A_vv);
ids = find(rows < cols);
edges = [rows(ids), cols(ids)]';
numE = size(edges, 2);
%
cols_A_vf = ones(3,1)*(1:numF);
rows_A_vf = Mesh.faceVIds;
vals_A_vf = ones(3, numF);
A_vf = sparse(rows_A_vf, cols_A_vf, vals_A_vf, numV, numF);
for vId = 1 : numV
    nFIds{vId} = find(A_vf(vId,:));
end
%
adjFaces = zeros(2, numE);
for eId = 1 : numE
    v1Id = edges(1, eId);
    v2Id = edges(2, eId);
    stats = [];
    for j = 1 : length(nFIds{v1Id})
        fId = nFIds{v1Id}(j);
        flag = is_face(Mesh.faceVIds(:, fId), [v1Id, v2Id]);
        if flag == 1
            stats = [stats, fId];
        end
    end
    if length(stats) ~= 2 && length(stats) ~= 1
        length(stats)
    end
    adjFaces(1:length(stats), eId) = stats;
end
%
tp1 = ones(2,1)*(1:numE);
tp2 = adjFaces;
tp3 = ones(2, numE);
tp1 = reshape(tp1, [1, 2*numE]);
tp2 = reshape(tp2, [1, 2*numE]);
tp3 = reshape(tp3, [1, 2*numE]);
ids = find(tp2 > 0);
tp1 = tp1(ids);
tp2 = tp2(ids);
tp3 = tp3(ids);
TP = sparse(tp1, tp2, tp3, numE, numF);
faceEIds = zeros(3, numF);
for fId = 1 : numF
    tp = find(TP(:,fId));
    reverse_ids = find(adjFaces(2, tp) == fId);
    faceEIds(:, fId) = tp';
    faceEIds(reverse_ids, fId) = -faceEIds(reverse_ids, fId);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag] = is_face(vids, edge)
flag = 0;
if min(vids(1:2)) == edge(1) && max(vids(1:2)) == edge(2)
    flag = 1;
    return;
end
if min(vids([1,3])) == edge(1) && max(vids([1,3])) == edge(2)
    flag = 1;
    return;
end
if min(vids([2,3])) == edge(1) && max(vids([2,3])) == edge(2)
    flag = 1;
    return;
end