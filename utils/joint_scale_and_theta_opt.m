function [theta_opt, scale_opt] = joint_scale_and_theta_opt(Mesh, theta_init, scale_init, lambda)
%
theta_opt = theta_init;
scale_opt = scale_init;
%
faceVPos1 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:));
faceVPos2 = Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
faceVPos3 = Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
edges1 = [Mesh.faceVIds(1,:);Mesh.faceVIds(2,:)];
edges2 = [Mesh.faceVIds(2,:);Mesh.faceVIds(3,:)];
edges3 = [Mesh.faceVIds(3,:);Mesh.faceVIds(1,:)];
%
edgeDiff1 = sum(Mesh.faceVFs.*(faceVPos1 - faceVPos2));
edgeDiff2 = sum(Mesh.faceVFs.*(faceVPos2 - faceVPos3));
edgeDiff3 = sum(Mesh.faceVFs.*(faceVPos3 - faceVPos1));
%
edges = [edges1, edges2, edges3];
numV = size(Mesh.vertexPoss, 2);
numE = size(edges, 2);
numF = size(Mesh.faceVIds, 2);
%
scaleIds = [1:numF,1:numF,1:numF];
edgeDiff_unscaled = [edgeDiff1, edgeDiff2, edgeDiff3];
% Regularization term
L_reg = face_laplacian(Mesh);
%
scale_opt = scale_init;
theta_opt = theta_init;
%
for iter = 1 : 64
    edgeDiff = edgeDiff_unscaled.*scale_opt(scaleIds);
    theta1 = theta_opt(edges(1,:));
    theta2 = theta_opt(edges(2,:)) + edgeDiff;
    vec_a = cos(theta1);
    vec_b = sin(theta1);
    vec_c = cos(theta2);
    vec_d = sin(theta2);
    diff1 = vec_a - vec_c;
    diff2 = vec_b - vec_d;
    g = double([diff1,diff2]');
    valsJ = [[-vec_b;vec_d;vec_d.*edgeDiff_unscaled],[vec_a;-vec_c;-vec_c.*edgeDiff_unscaled]];
    colsJ = [[edges;scaleIds+numV],[edges;scaleIds+numV]];
    rowsJ = ones(3,1)*(1:(2*numE));
    J = sparse(rowsJ, colsJ, double(valsJ), 2*numE, numV + numF);
    A = J'*J;
    % Update the regularization term
    ids = (numV+1):(numV+numF);
    A(ids,ids) = A(ids,ids) + lambda*L_reg;
    b =-J'*g;
    v = sparse([ones(1,numV), zeros(1,numF)]);
    % Augmented Lagragian
    A = [A,v';v, 0];
    b = [b;0];
    dx_cur = A\b;
    dx_cur = dx_cur(1:(numV+numF));
    x_cur = [theta_opt, scale_opt]';
    e_cur = objective_value(edges, numV, x_cur, scaleIds, edgeDiff_unscaled, L_reg, lambda);
    alpha = 1;
    for searchId = 1 : 20
 %       if min(x_cur(ids)+alpha*dx_cur(ids)) < 0
 %           alpha = alpha/2;
 %           continue;
 %       end
        e_next = objective_value(edges, numV, x_cur+alpha*dx_cur, scaleIds, edgeDiff_unscaled, L_reg, lambda);
        if e_next < e_cur
            theta_opt = theta_opt+alpha*dx_cur(1:numV)';
            scale_opt = scale_opt+alpha*dx_cur(ids)';
            break;
        end
        alpha = alpha/2;
    end
    fprintf('e_next = %f.\n', e_next);
end
%
function [e] = objective_value(edges, numV, x_cur, scaleIds, edgeDiff_unscaled, L_reg, lambda)
%
theta = x_cur(1:numV)';
scale = x_cur((numV+1):length(x_cur))';
dif = theta(edges(1,:)) - theta(edges(2,:)) - edgeDiff_unscaled.*scale(scaleIds);
dcos = cos(dif) - 1;
dsin = sin(dif);
e = sum(dcos.*dcos) + sum(dsin.*dsin) + lambda*(scale*L_reg*scale');
%
function [L_reg] = face_laplacian(Mesh)
%
numF = size(Mesh.faceVIds, 2);
tp = find(min(Mesh.adjFaces) > 0);
tp = Mesh.adjFaces(:, tp);
L_reg = sparse(tp(1,:), tp(2,:), ones(1,length(tp)), numF, numF);
L_reg = L_reg + L_reg';
d = full(sum(L_reg));
L_reg = sparse(1:length(d),1:length(d),d) - L_reg;