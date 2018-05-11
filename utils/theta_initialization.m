function [theta_init] = theta_initialization(Mesh, scale_fixed)
% We use spectral initialization to obtain the initial theta
faceVPos1 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:));
faceVPos2 = Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
faceVPos3 = Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
faceVFs = Mesh.faceVFs.*(ones(3,1)*scale_fixed);
%
edges1 = [Mesh.faceVIds(1,:);Mesh.faceVIds(2,:)];
edges2 = [Mesh.faceVIds(2,:);Mesh.faceVIds(3,:)];
edges3 = [Mesh.faceVIds(3,:);Mesh.faceVIds(1,:)];
%
edgeDiff1 = sum(faceVFs.*(faceVPos1 - faceVPos2));
edgeDiff2 = sum(faceVFs.*(faceVPos2 - faceVPos3));
edgeDiff3 = sum(faceVFs.*(faceVPos3 - faceVPos1));
%
edges = [edges1, edges2, edges3];
edgeDiff = [edgeDiff1, edgeDiff2, edgeDiff3];
numE = size(edges, 2);
VALs = zeros(2, 2*numE);
%
VALs(1, 1:2:(2*numE)) = cos(edgeDiff);
VALs(1, 2:2:(2*numE)) = -sin(edgeDiff);
VALs(2, 1:2:(2*numE)) = sin(edgeDiff);
VALs(2, 2:2:(2*numE)) = cos(edgeDiff);
%
ROWs = [kron(2*edges(1,:)-1, ones(1,2));kron(2*edges(1,:), ones(1,2))];
%
COLs(:, 1:2:(2*numE)) = ones(2,1)*(2*edges(2,:)-1);
COLs(:, 2:2:(2*numE)) = ones(2,1)*(2*edges(2,:));
numV = size(Mesh.vertexPoss, 2);
Adj_connection = sparse(ROWs, COLs, VALs, 2*numV, 2*numV);
Adj_connection = Adj_connection + Adj_connection';
%
Adj = sparse(edges(1,:), edges(2,:), ones(1,size(edges,2)), numV, numV);
Adj = Adj + Adj';
d = sum(full(Adj));
Lap_connection = sparse(1:(2*numV),1:(2*numV), kron(d,ones(1,2))) - Adj_connection;
[eigenVecs, eigenVals] = eigs(Lap_connection, 2, 1e-10);
angles = reshape(eigenVecs(:,1), [2, numV]);
norms = sqrt(sum(angles.*angles));
angles = angles./(ones(2,1)*norms);
theta_init = acos(angles(1,:));
ids = find(angles(2,:) < 0);
theta_init(ids) = - theta_init(ids);