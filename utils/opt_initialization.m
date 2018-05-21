function [scale_init, func_init] = opt_initialization(Mesh, basis_scale)
% J*x - W*s == dif
OPT = mesh_opt_structure(Mesh);
scale_init = basis_scale.eigVecs(:,1);
if sum(scale_init) < 0
    scale_init = -scale_init;
end
dif = OPT.W*scale_init;
[dim, numV] = size(OPT.J);
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

function [OPT] = mesh_opt_structure(Mesh)
%
numF = size(Mesh.faceEIds, 2);
numE = size(Mesh.edges, 2);
numV = size(Mesh.vertexPoss,2);
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
c1 = scales_1;
fids_1 = Mesh.adjFaces(1,ids1)';
c2 = scales_2;
fids_2 = Mesh.adjFaces(2,ids2)';
OPT.J = [J1; J2];
OPT.c = [c1, c2];
OPT.fids = [fids_1; fids_2];
dim = size(OPT.J, 1);
OPT.W = sparse(1:dim, OPT.fids, double(OPT.c), dim, numF);

function [faceNormals] = mesh_face_normal(Mesh)
%
e12 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:)) - Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
e13 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:)) - Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
faceNormals = cross(e12, e13);
norms = sqrt(sum(faceNormals.*faceNormals));
faceNormals(1,:) = faceNormals(1,:)./norms;
faceNormals(2,:) = faceNormals(2,:)./norms;
faceNormals(3,:) = faceNormals(3,:)./norms;

