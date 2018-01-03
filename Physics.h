#ifndef PHYSICS_H
#define PHYSICS_H

#include <Eigen/Core>
#include <vector>

struct MeshData;

// precomputed data about the mesh needed for the physics
struct PhysicsData
{
    std::vector<Eigen::Matrix2d> as; // first fundamental form on each face
    std::vector<Eigen::Matrix2d> bs; // second fundamental form on each face
    std::vector<Eigen::MatrixXd> das; // 3x9 matrix of derivatives; each row maps from changes in a face's vertices to changes in an entry of the first fundamental form (unrolled as [a_11 a_12=a_21 a_22]^T)
    std::vector<Eigen::MatrixXd> dbs; // 3x18 matrix of derivatives; each row maps from changes in a face's vertices and the "spoke" vertex of the neighboring triangles, 
                                      // to changes in an entry of the second fundamental form (unrolled as [b_11 b_12=b_21 b_22]^T)
};

void physicsDataFromMesh(const MeshData &mesh, PhysicsData &phys);

void physicalForces(const MeshData &mesh, const PhysicsData &fundamentalForms, const Eigen::VectorXd &widthFractions, const Eigen::VectorXd &vectorField, double kb, double kt,
    Eigen::MatrixXd &forceField);

void testForces(const MeshData &mesh, const PhysicsData &fundamentalForms, const Eigen::VectorXd &widthFractions, const Eigen::VectorXd &vectorField, double kb, double kt, int vertex, int coord);

#endif
