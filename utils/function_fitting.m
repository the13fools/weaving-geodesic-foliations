function [scale_init, func_init] = function_fitting(Mesh, lambda)
% J*x - W*s == dif
numF = size(Mesh.faceEIds, 2);
numE = size(Mesh.edges, 2);
numV = size(Mesh.vertexPoss, 2);
valsJ_fe = 2*(Mesh.faceEIds > 0) - 1;
rowsJ_fe = ones(3,1)*(1:numF);
colsJ_fe = abs(Mesh.faceEIds);
OPT.J_fe = sparse(rowsJ_fe, colsJ_fe, valsJ_fe, numF, numE);
%
ids1 = find(Mesh.adjFaces(1,:) > 0);
ids2 = find(Mesh.adjFaces(2,:) > 0);
J_ev = sparse(ones(2,1)*(1:numE), Mesh.edges, [1,-1]'*ones(1, numE));
edgeVec = Mesh.vertexPoss(:, Mesh.edges(1,:))...
    - Mesh.vertexPoss(:, Mesh.edges(2,:));
scales_1 = sum(edgeVec(:, ids1).*Mesh.faceVFs(:, Mesh.adjFaces(1,ids1)));
scales_2 = sum(edgeVec(:, ids2).*Mesh.faceVFs(:, Mesh.adjFaces(2,ids2)));
%
J1 = J_ev(ids1,:);
J2 = J_ev(ids2,:);
eids_1 = ids1;
eids_2 = ids2;
c1 = scales_1;
fids_1 = Mesh.adjFaces(1,ids1)';
c2 = scales_2;
fids_2 = Mesh.adjFaces(2,ids2)';
OPT.J = [J1; J2];
OPT.c = [c1, c2];
OPT.fids = [fids_1; fids_2];
dim = size(OPT.J, 1);
OPT.W = sparse(1:dim, OPT.fids, double(OPT.c), dim, numF);
L = OPT.W'*OPT.W - (OPT.J'*OPT.W)'*pinv(full(OPT.J'*OPT.J))*(OPT.J'*OPT.W);
[i,j,s] = find(OPT.W'*OPT.W);
% The Laplacian on the mesh faces
L_reg = face_laplacian(Mesh);
L = L + lambda*L_reg;
[u,v] = eigs(sparse(L), 1, 1e-10);
s = u;
if sum(s) < 0
    s = -s;
end
scale_init = s/mean(s);
dif = OPT.W*scale_init;
%
cvx_begin
variable func_init(numV);
variable y(dim);
minimize (y'*y)
subject to
y == OPT.J*func_init - dif;
sum(func_init) == 0;
cvx_end
scale_init = scale_init';
func_init = func_init';

function [L_reg] = face_laplacian(Mesh)
%
numF = size(Mesh.faceVIds, 2);
tp = find(min(Mesh.adjFaces) > 0);
tp = Mesh.adjFaces(:, tp);
L_reg = sparse(tp(1,:), tp(2,:), ones(1,length(tp)), numF, numF);
L_reg = L_reg + L_reg';
d = full(sum(L_reg));
L_reg = sparse(1:length(d),1:length(d),d) - L_reg;

function [faceNormals] = mesh_face_normal(Mesh)
%
e12 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:)) - Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
e13 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:)) - Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
faceNormals = cross(e12, e13);
norms = sqrt(sum(faceNormals.*faceNormals));
faceNormals(1,:) = faceNormals(1,:)./norms;
faceNormals(2,:) = faceNormals(2,:)./norms;
faceNormals(3,:) = faceNormals(3,:)./norms;

