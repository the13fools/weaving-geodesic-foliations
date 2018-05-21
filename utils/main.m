% parameters
lambda_initializing_scale = 1e1;
lambda_almin_scale = 1e-2;
% Paths of file names
cut_mesh_file = 'debug_cut.obj';
ori_mesh_file = 'debug_augmented.obj';
vector_field_file = 'debug.field';
% load meshes
cut_mesh = load_obj(cut_mesh_file);
ori_mesh = load_obj(ori_mesh_file);
% computer mesh topology
[cut_mesh.edges, cut_mesh.adjFaces, cut_mesh.faceEIds] = mesh_topology(cut_mesh);
[ori_mesh.edges, ori_mesh.adjFaces, ori_mesh.faceEIds] = mesh_topology(ori_mesh);
% load vector field
unoriented_vfs = load(vector_field_file)';
% rotate the vector field
facenormals = mesh_face_normal(cut_mesh);
cut_mesh.faceVFs = cross(facenormals, unoriented_vfs);
ori_mesh.faceVFs = cut_mesh.faceVFs;
%
basis_scale_cut_mesh = reduced_basis_scale(cut_mesh, lambda_initializing_scale, 1);
%
[scale_init, func_init] = opt_initialization(ori_mesh, basis_scale_cut_mesh);
% A heuristic for scaling scale_init and func_init
[maxEdgeDiff, edgeId] = max(abs(func_init(ori_mesh.edges(1,:))-func_init(ori_mesh.edges(2,:))));
scale = .5*pi/maxEdgeDiff;
% Apply the scaling
scale_opt = scale_init*scale;
func_opt = func_init*scale;
% Initialize theta
theta_opt = theta_initialization(ori_mesh, scale_opt);
% Alternating minimization
for iter = 1:128
    scale_opt = optimize_scales(ori_mesh, theta_opt, scale_opt, lambda_almin_scale);
    theta_opt = optimize_function(ori_mesh, theta_opt, scale_opt);
end

trisurf(ori_mesh.faceVIds', ori_mesh.vertexPoss(1,:), ori_mesh.vertexPoss(2,:), ori_mesh.vertexPoss(3,:), cos(theta_opt))
daspect([1 1 1])

%%%%%% Visualize

numfaces = size(ori_mesh.faceVIds, 2)
numverts = size(ori_mesh.vertexPoss, 2)

shift = repmat(1:6,[numverts/6 1]);
shift = shift(:)';

ori_mesh.vertexPossShift = [ ori_mesh.vertexPoss(1,:) + shift; ori_mesh.vertexPoss(2,:); ori_mesh.vertexPoss(3,:) ];
ori_mesh.shiftEdgeLens = [ sum( ( ori_mesh.vertexPossShift(:, ori_mesh.faceVIds(1, :)) - ori_mesh.vertexPossShift(:, ori_mesh.faceVIds(2, :)) ).^2, 1 );
    sum( ( ori_mesh.vertexPossShift(:, ori_mesh.faceVIds(2, :)) - ori_mesh.vertexPossShift(:, ori_mesh.faceVIds(3, :)) ).^2, 1 );
    sum( ( ori_mesh.vertexPossShift(:, ori_mesh.faceVIds(3, :)) - ori_mesh.vertexPossShift(:, ori_mesh.faceVIds(1, :)) ).^2, 1 ); ];
ori_mesh.faceLongEdges = max(ori_mesh.shiftEdgeLens, [],1);

trisurf(ori_mesh.faceVIds(:, find(ori_mesh.faceLongEdges < .1 ))', ori_mesh.vertexPossShift(1,:), ori_mesh.vertexPossShift(2,:), ori_mesh.vertexPossShift(3,:), cos(theta_opt))

trisurf(ori_mesh.faceVIds', Torus2_ori.vertexPoss(1,:) + shift, Torus2_ori.vertexPoss(2,:), Torus2_ori.vertexPoss(3,:), cos(theta_altermin_opt))



%%%%%%  Save

scale_opt_save = scale_opt';
theta_opt_save = theta_opt';
save("anuerism_small.sval", "scale_opt_save", '-ascii');
save("anuerism_small.theta", "theta_opt_save", '-ascii');


