function [basis_scale] = reduced_basis_scale(Mesh, lambda, dim_basis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J*x - W*s == dif
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

% OPT.X = rand(size(OPT.W, 2), 1);
% 
OPT.lap = OPT.J'*OPT.J;
OPT.lap2 = [OPT.lap, ones(numV,1);
    ones(1,numV), 0];
OPT.WtW = OPT.W'*OPT.W;
OPT.JtW = OPT.J'*OPT.W;
%
L_reg = face_laplacian(Mesh);
maxL = max(full(sum(abs(OPT.WtW)))) + lambda*max(full(sum(abs(L_reg))))/2;
% OPT.afun = @(x) OPT.W'*x*OPT.W - OPT.W'*OPT.J * ( OPT.lap \( OPT.J' * x ) )*OPT.W + lambda * L_reg*x;
OPT.afun = @(x) maxL*x - (OPT.WtW*x - OPT.JtW'*bfunc(OPT.lap2, OPT.JtW*x) + lambda * L_reg *x);
[basis_scale.eigVecs, basis_scale.eigVals] = eigs(OPT.afun, numF, dim_basis, 'LM');
basis_scale.eigVals = maxL - diag(basis_scale.eigVals)';
if 0
    A = full(OPT.WtW - OPT.JtW'*pinv(full(OPT.lap))*OPT.JtW + lambda*L_reg);
    [u,v] = eig((A+A')/2);
    v = diag(v)';
end
function [lap2x] = bfunc(Lap2, x)
y = Lap2\[x;0];
lap2x = y(1:(length(y)-1));

function [L_reg] = face_laplacian(Mesh)
%
numF = size(Mesh.faceVIds, 2);
tp = find(min(Mesh.adjFaces) > 0);
tp = Mesh.adjFaces(:, tp);
L_reg = sparse(tp(1,:), tp(2,:), ones(1,length(tp)), numF, numF);
L_reg = L_reg + L_reg';
d = full(sum(L_reg));
L_reg = sparse(1:length(d),1:length(d),d) - L_reg;


