function [scale_opt] = optimize_scales(Mesh, theta_fixed, scale_init, lambda)
% Solve for the initial theta
%
numF = size(Mesh.faceVIds, 2);
faceVPos1 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:));
faceVPos2 = Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
faceVPos3 = Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
scaleIds = [1:numF,1:numF,1:numF];
edges1 = [Mesh.faceVIds(1,:);Mesh.faceVIds(2,:)];
edges2 = [Mesh.faceVIds(2,:);Mesh.faceVIds(3,:)];
edges3 = [Mesh.faceVIds(3,:);Mesh.faceVIds(1,:)];
%
edgeDiff1 = sum(Mesh.faceVFs.*(faceVPos1 - faceVPos2));
edgeDiff2 = sum(Mesh.faceVFs.*(faceVPos2 - faceVPos3));
edgeDiff3 = sum(Mesh.faceVFs.*(faceVPos3 - faceVPos1));
%
edges = [edges1, edges2, edges3];
numE = size(edges, 2);
numV = size(Mesh.vertexPoss, 2);
edgeDiff_unscaled = [edgeDiff1, edgeDiff2, edgeDiff3];
targetDiff = theta_fixed(edges(1,:)) - theta_fixed(edges(2,:));
vec_a = cos(targetDiff);
vec_b = sin(targetDiff);
scale_opt = scale_init;
% Regularization term
L_reg = face_laplacian(Mesh);
for iter = 1 : 8
    edgeDiff = edgeDiff_unscaled.*scale_opt(scaleIds);
    vec_c = cos(edgeDiff);
    vec_d = sin(edgeDiff);
    diff1 = vec_a - vec_c;
    diff2 = vec_b - vec_d;
    g1 = [diff1,diff2]';
    valsJ1 = [vec_d.*edgeDiff_unscaled,-vec_c.*edgeDiff_unscaled]';
    colsJ1 = [scaleIds,scaleIds]';
    J1 = sparse(1:length(colsJ1), colsJ1, double(valsJ1));
    A = J1'*J1 + lambda*L_reg;
    b =-J1'*double(g1);
    dscale = A\b;
    e_cur = objective_value(vec_a, vec_b, L_reg, lambda, scale_opt, edgeDiff_unscaled, scaleIds);
    alpha = 1;
    for searchId = 1 : 20
        e_next = objective_value(vec_a, vec_b, L_reg, lambda, scale_opt+alpha*dscale',...
            edgeDiff_unscaled, scaleIds);
        if e_next < e_cur
            scale_opt = scale_opt+alpha*dscale';
            break;
        end
        alpha = alpha/2;
    end
end
fprintf('e_next = %f\n', e_cur);

function [e] = objective_value(vec_a, vec_b, L_reg, lambda, scale_cur, edgeDiff_unscaled, scaleIds)
%
edgeDiff = edgeDiff_unscaled.*scale_cur(scaleIds);
vec_c = cos(edgeDiff);
vec_d = sin(edgeDiff);
diff1 = vec_a - vec_c;
diff2 = vec_b - vec_d;
e = diff1*diff1' + diff2*diff2' + lambda*(scale_cur*L_reg*scale_cur');
    
function [L_reg] = face_laplacian(Mesh)
%
numF = size(Mesh.faceVIds, 2);
tp = find(min(Mesh.adjFaces) > 0);
tp = Mesh.adjFaces(:, tp);
L_reg = sparse(tp(1,:), tp(2,:), ones(1,length(tp)), numF, numF);
L_reg = L_reg + L_reg';
d = full(sum(L_reg));
L_reg = sparse(1:length(d),1:length(d),d) - L_reg;
