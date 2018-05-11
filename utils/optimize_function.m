function [theta_opt] = optimize_function(Mesh, theta_init, scale_fixed)
% Solve for the initial theta
theta_opt = theta_init - mean(theta_init);
%
faceVPos1 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:));
faceVPos2 = Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
faceVPos3 = Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
faceVFs = Mesh.faceVFs.*(ones(3,1)*scale_fixed);
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
numV = size(Mesh.vertexPoss, 2);
for iter = 1 : 8
    vals_left = theta_opt(edges(1,:));
    vals_right = theta_opt(edges(2,:)) + edgeDiff;
    vec_a = cos(vals_left);
    vec_b = sin(vals_left);
    vec_c = cos(vals_right);
    vec_d = sin(vals_right);
    g = reshape([vec_a - vec_c;vec_b - vec_d], [2*numE,1]);
    valsJ = zeros(2*numE, 2);
    valsJ(1:2:(2*numE),1) = -vec_b';
    valsJ(2:2:(2*numE),1) = vec_a';
    valsJ(1:2:(2*numE),2) = vec_d';
    valsJ(2:2:(2*numE),2) = -vec_c';
    rowsJ = (1:(2*numE))'*ones(1,2);
    colsJ = zeros(2*numE, 2);
    colsJ(1:2:(2*numE),1) = edges(1,:)';
    colsJ(2:2:(2*numE),1) = edges(1,:)';
    colsJ(1:2:(2*numE),2) = edges(2,:)';
    colsJ(2:2:(2*numE),2) = edges(2,:)';
    J = sparse(rowsJ, colsJ,valsJ, 2*numE, numV);
    A = J'*J;
    b =-J'*double(g);
    H = [A, ones(numV,1);ones(1,numV),0];
    h = [b; 0];
    dtheta = H\h;
    dtheta = dtheta(1:numV);
    e_cur = objective_value(edges, edgeDiff, theta_opt);
    alpha = 1;
    for searchId = 1 : 10
        e_next = objective_value(edges, edgeDiff, theta_opt+alpha*dtheta');
        if e_next < e_cur
            theta_opt = theta_opt+alpha*dtheta';
            break;
        end
        alpha = alpha/2;
    end
end

function [e] = objective_value(edges, edgeDiff, theta_opt)
dif = theta_opt(edges(1,:)) - edgeDiff - theta_opt(edges(2,:));
difx = cos(dif) - 1;
dify = sin(dif);
e = difx*difx' + dify*dify';
