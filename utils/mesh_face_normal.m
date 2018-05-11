function [faceNormals] = mesh_face_normal(Mesh)
%
e12 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:)) - Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
e13 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:)) - Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
faceNormals = cross(e12, e13);
norms = sqrt(sum(faceNormals.*faceNormals));
faceNormals(1,:) = faceNormals(1,:)./norms;
faceNormals(2,:) = faceNormals(2,:)./norms;
faceNormals(3,:) = faceNormals(3,:)./norms;