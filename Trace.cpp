#include "Trace.h"
#include "VectorUtils.h"
#include "Distance.h"
#include "Weave.h"

#include <fstream>
#include <iostream>

// For making directories
#include <sys/stat.h>
 // #include <direct.h>

Collision::Collision(int rod1, int rod2, int seg1, int seg2) : rod1(rod1), rod2(rod2), 
                                                               seg1(seg1), seg2(seg2)
{
} 

Trace::Trace(){}
Trace::~Trace(){}

void Trace::popLastCurve()
{
    if (curves.size() > 0) 
    {
        curves.pop_back();
	normals.pop_back();
	modes.pop_back();
        bending.pop_back();
    }
}

Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2)
{
    Eigen::Vector3d t1 = e1 / e1.norm();
    Eigen::Vector3d t2 = e2 / e2.norm();
    Eigen::Vector3d n = t1.cross(t2);
    if (n.norm() < 1e-8)
        return v;
    n /= n.norm();
    Eigen::Vector3d p1 = n.cross(t1);
    Eigen::Vector3d p2 = n.cross(t2);
    return v.dot(n)*n + v.dot(t1)*t2 + v.dot(p1)*p2;
}


void Trace::loadSampledCurves(const std::string &filename)
{
    std::ifstream curve_file("f_curve_xyz.txt");
    if (!curve_file)
    {
        std::cerr << "Couldn't load trace file" << std::endl;
        return;
    }
    std::ifstream normal_file("face_normal.txt");
    if (!normal_file)
    {
        std::cerr << "Couldn't load trace file" << std::endl;
        return;
    }

    int npoints, dum1, dum2;
    double dum3, dum4, dum5;
    while (normal_file >> dum3 >> dum4 >> dum5)
    {
        curve_file >> npoints >> dum1 >> dum2;
        Eigen::MatrixXd curve = Eigen::MatrixXd::Zero(npoints, 3);
        Eigen::MatrixXd normal = Eigen::MatrixXd::Zero(npoints, 3);
        Eigen::VectorXd bend = Eigen::VectorXd::Zero(npoints); // Not implemented yet

        Eigen::Vector3d next_normal(0, 0, 0);
        for (int i = 0; i < npoints; i++)
        {
            curve_file >> curve(i, 0) >> curve(i, 1) >> curve(i, 2);
            normal_file >> normal(i, 0) >> normal(i, 1) >> normal(i, 2);
            Eigen::Vector3d tmp = normal.row(i);
            if (tmp.dot(next_normal) < 0.)
            {
                tmp = -tmp;
            }
            normal.row(i) = next_normal;
            next_normal = tmp;
        }

        std::cout << normal.row(npoints - 1) << std::endl;



        if (npoints > 100)
        {
            curves.push_back(curve);
            normals.push_back(normal);
            bending.push_back(bend);
            modes.push_back(Trace_Mode::FIELD); // Make another render option
        }
    }
    curve_file.close();
    normal_file.close();
}



void Trace::logRibbonsToFile(std::string foldername, std::string filename)
{
    std::stringstream folderpath;
    folderpath << "log/" << foldername;
#ifndef WIN32
    mkdir("log", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(folderpath.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
    folderpath << filename << ".rod";
    std::ofstream myfile(folderpath.str());

    if (!myfile.is_open())
    {
        std::cout << "Unable to open file";
        return;
    }

    // Find collisions between rods
    std::vector<Collision> collisions;
    for (int i = 0; i < curves.size(); i++)
    {
        for (int j = i; j < curves.size(); j++)
        {
            computeIntersections(i, j, collisions);
        }
    }

    // Write Header 
    myfile << curves.size() << std::endl;;
    myfile << collisions.size() << std::endl;;
  //  myfile << 0 << std::endl;;
    myfile << "0.001"  << std::endl;;
    myfile << "1e+08"  << std::endl;;
    myfile << "1"  << std::endl  << std::endl  << std::endl;;


    // Decimate and log rods
    std::vector<Eigen::VectorXd> desc_maps;
    std::vector< std::vector<Eigen::Vector3d> > desc_curves; //eeew
    for (int curveId = 0; curveId < curves.size(); curveId++)
    {
        Eigen::MatrixXd curve = curves[curveId];
        Eigen::MatrixXd curveNormals = normals[curveId];

        double max_length = 0.;
        for (int i = 0; i < curve.rows() - 1; i++)
        {
            double seg_length = (curve.row(i) - curve.row(i + 1)).norm();
            if (seg_length > max_length)
            {
                max_length = seg_length;
            }
        }

        // Decimate 
        std::vector<Eigen::Vector3d> cnew;
        std::vector<Eigen::Vector3d> nnew;
        std::vector<double> widths;
        cnew.push_back(curve.row(0));
        int seg_counter = 0;
        Eigen::VectorXd desc_mapping = Eigen::VectorXd::Zero(curve.rows());
        Eigen::Vector3d prev_point = cnew.back();
        for (int i = 1; i < curve.rows(); i++)
        {
            Eigen::Vector3d curr_point = curve.row(i);
            double seg_length = (prev_point - curr_point).norm();
            desc_mapping(i-1) = seg_counter;
            if (seg_length > max_length)
            {
                seg_counter++;
                cnew.push_back(curve.row(i));
                Eigen::Vector3d currEdge = curve.row(i) - curve.row(i - 1);
                Eigen::Vector3d targEdge = cnew[seg_counter] - cnew[seg_counter - 1];
                nnew.push_back(parallelTransport(curveNormals.row(i), currEdge, targEdge));
                widths.push_back(.005);
                prev_point = cnew.back();
            }
        }
//        cnew.pop_back();
//        cnew.push_back(curve.row(curve.rows() - 1));


        desc_maps.push_back(desc_mapping);
        desc_curves.push_back(cnew);


/*    Eigen::Vector3d prev(0,0,0); 
    Eigen::Vector3d cur(0,0,0); 
    for (int i = 1; i < cnew.size(); i++) 
    {
	Eigen::Vector3d n = nnew[i-1];
	prev = cnew[i-1];
	cur =  cnew[i];
        std::cout << (cur-prev).dot(n) << std::endl;
    }
*/
        if (cnew.size() < 3) 
            continue;


std::cout << "blah";
        myfile << cnew.size() << "\n";
        myfile << "0\n";
        for (int i = 0; i < cnew.size(); i++)
        {
            myfile << cnew[i](0) << " " << cnew[i](1) << " " << cnew[i](2) << " ";
        }
        myfile << "\n";

std::cout << "blah1" << std::endl;
//        myfile << nnew[0](0) << " " << nnew[0](1) << " " << nnew[0](2) << " ";
//std::cout << "blah2";
        for (int i = 0; i < nnew.size(); i++)
        {
            myfile << nnew[i](0) << " " << nnew[i](1) << " " << nnew[i](2) << " ";
        }
        myfile << "\n";

        for (int i = 0; i < widths.size(); i++)
        {
            myfile << " " << widths[i];
        }
        myfile << "\n";

    }
    for (int i = 0; i < collisions.size(); i++)
    {
        Collision col = collisions[i];
        std::vector<Eigen::Vector3d> c1 = desc_curves[col.rod1];
        std::vector<Eigen::Vector3d> c2 = desc_curves[col.rod2];
        int idx1 = desc_maps[col.rod1](col.seg1);
        int idx2 = desc_maps[col.rod2](col.seg2);

        //	std::cout << "r1 " << col.rod1 << "r2 " << col.rod2 << "idx1 " << idx1 << "idx2 " << idx2 <<"\n";
        double p0bary, p1bary, q0bary, q1bary;

        Eigen::Vector3d dist = Distance::edgeEdgeDistance(c1[idx1],
            c1[idx1 + 1],
            c2[idx2],
            c2[idx2 + 1],
            p0bary, p1bary, q0bary, q1bary);

        //	std::cout << "dist " << dist.norm() << "idx1" << idx1 << "idx2" << idx2;
        myfile << col.rod1 << " " << col.rod2 << " " << idx1 << " " << idx2
            << " " << p1bary << " " << q1bary << " " << 1000. << "\n";

    }


    myfile.close();
}


// This: Recieves a point on an edge of a triangle
//       Returns the vertex id opposing that edge in the triangle
int getOpVId(const Weave &wv, const Eigen::Vector3d prev_point, int faceId)
{
    for (int j = 0; j < 3; j++)
    {
        if(wv.faceNeighbors(faceId, j) == -1)
            continue;
	Eigen::VectorXi e = wv.edgeVerts.row(wv.faceNeighbors(faceId, j));
	Eigen::Vector3d e_test =  wv.V.row( e(0) ) - wv.V.row( e(1) );
	e_test.normalize();
        Eigen::Vector3d point_test = wv.V.row( e(0) );
        point_test = ( prev_point - point_test );
        point_test.normalize();
	if ( e_test.dot( point_test ) > .99 )
	{
	    Eigen::VectorXi e_next = wv.edgeVerts.row(wv.faceNeighbors(faceId, (j + 1) % 3 ));
	    if (e_next(0) == e(0) || e_next(0) == e(1) )
	    {
		    return e_next(1);
	    }
	    else 
	    {
		    return e_next(0);
	    }
	} 
    }
    return -1;
}

int getOpVIdFromEdge(const Weave &wv, int curr_edge, int faceId)
{
    Eigen::Vector2i e = wv.edgeVerts.row(curr_edge);
    for (int i = 0; i < 3; i++)
    {
        if (wv.F(faceId, i) != e(0) && wv.F(faceId, i) != e(1))
        {
            return wv.F(faceId, i);
        }
    }
    assert(false);
    return -1;
}

void Trace::computeIntersections(int curveIdx1, int curveIdx2, std::vector<Collision> &collisions)
{
    Eigen::MatrixXd c1 = curves[curveIdx1];
    Eigen::MatrixXd c2 = curves[curveIdx2];

    for (int i = 0; i < c1.rows() - 1; i++)
    {
        for (int j = 0; j < c2.rows() - 1; j++)
        {
            if (curveIdx1 == curveIdx2 && (i - j) < 2) { continue; }
            else
            {
                double p0bary, p1bary, q0bary, q1bary;
                Eigen::Vector3d dist = Distance::edgeEdgeDistance(c1.row(i),
                    c1.row(i + 1),
                    c2.row(j),
                    c2.row(j + 1),
                    p0bary, p1bary, q0bary, q1bary);
                if (dist.norm() < 1e-6 && p0bary != 0 && p0bary != 1.0 && q0bary != 0 && q0bary != 1.0)
                {
                    collisions.push_back(Collision(curveIdx1, curveIdx2, i, j));
                }
            }
        }
    }
}

/*
void Trace::projectToClosestFaceDir()
{

}
*/
/*
Eigen::VectorXd mapFaceVectorToGlobalCoordinates(const Weave &wv, const Eigen::VectorXi coords, int faceId)
{
    Eigen::VectorXd ret = Eigen::VectorXd::Zero(3);
    for (int i = 0; i < 3; i++)
    {
        
    }
}
*/

void Trace::traceCurve(const Weave &wv, const Trace_Mode trace_state,
    int traceIdx, int sign, int faceId, int steps)
{
    Eigen::MatrixXd curve = Eigen::MatrixXd::Zero(steps, 3);
    Eigen::MatrixXd normal = Eigen::MatrixXd::Zero(steps, 3);
    curve.row(0) = (1. / 3. * wv.V.row(wv.F(faceId, 0)) +
        1. / 3. * wv.V.row(wv.F(faceId, 1)) +
        1. / 3. * wv.V.row(wv.F(faceId, 2)));
    normal.row(0) = (faceNormal(wv.F, wv.V, faceId));


    
    int curr_dir_idx = abs(traceIdx);
    double coeff_dir = sign > 0 ? 1 : -1;
    Eigen::Vector3d curr_dir = coeff_dir * wv.Bs[faceId] * wv.v(faceId, curr_dir_idx);

    // assumes roughly delaunay    
    curr_dir = curr_dir.normalized() * wv.averageEdgeLength * 1000.;

    int curr_face_id = faceId;

    // Project backwards to initialize.  Kinda hacky.
    double min_dist = std::numeric_limits<double>::infinity();
    Eigen::Vector3d startpoint;
    startpoint.setZero();
    for (int i = 0; i < 3; i++)
    {
        startpoint += wv.V.row(wv.F(faceId, i));
    }
    startpoint /= 3.0;

    int curr_edge_id = -1;
    Eigen::Vector3d prev_point;

    for (int i = 0; i < 3; i++)
    {
        int next_edge_id = wv.faceEdges(curr_face_id, i);
        Eigen::Vector3d op_v1 = wv.V.row(wv.edgeVerts(next_edge_id, 0));
        Eigen::Vector3d op_v2 = wv.V.row(wv.edgeVerts(next_edge_id, 1));

        double p0bary, p1bary, q0bary, q1bary;
        Eigen::Vector3d dist = Distance::edgeEdgeDistance(startpoint,
            startpoint - curr_dir,
            op_v1, op_v2,
            p0bary, p1bary, q0bary, q1bary);
        if (dist.norm() < min_dist)
        {
            min_dist = dist.norm();
            prev_point = op_v1 * q0bary + op_v2 * q1bary;
            curr_edge_id = next_edge_id;
        }
    }

    assert(curr_edge_id > -1);
    int op_v_id = getOpVIdFromEdge(wv, curr_edge_id, curr_face_id);
    Eigen::VectorXd bend = Eigen::VectorXd::Zero(steps);
    for (int i = 1; i < steps; i++)
    {
        Eigen::Vector3d op_vertex = wv.V.row(op_v_id);
        Eigen::Vector3d split = op_vertex - prev_point;

        double split_len = split.norm();
        split = split / split_len;
        Eigen::Vector3d n = faceNormal(wv.F, wv.V, curr_face_id);
        Eigen::Vector3d perp = split.cross(n);
        normal.row(i) = n;

        int op_edge_id = -1;
        for (int j = 0; j < 3; j++)
        {
            Eigen::Vector2i e = wv.edgeVerts.row(wv.faceEdges(curr_face_id, j));

            if (e(0) == op_v_id || e(1) == op_v_id)
            {
                Eigen::Vector3d e_test = (wv.V.row(e(0)) + wv.V.row(e(1))) * .5;
                e_test -= prev_point;
                if (e_test.dot(perp) * curr_dir.dot(perp) > 0.)
                {
                    op_edge_id = j;
                    break;
                }
            }
        }
        // stop if we hit vertex
        if (op_edge_id == -1)
        {
            curve.conservativeResize(i, 3);
            normal.conservativeResize(i, 3);
            break;
        }

        // Find intersection point.
        int next_edge_id = wv.faceEdges(curr_face_id, op_edge_id);
        Eigen::Vector3d op_v1 = wv.V.row(wv.edgeVerts(next_edge_id, 0));
        Eigen::Vector3d op_v2 = wv.V.row(wv.edgeVerts(next_edge_id, 1));

        double p0bary, p1bary, q0bary, q1bary;
        Eigen::Vector3d dist = Distance::edgeEdgeDistance(prev_point,
            prev_point + curr_dir,
            op_v1, op_v2,
            p0bary, p1bary, q0bary, q1bary);

        curve.row(i) = (op_v1 * q0bary + op_v2 * q1bary);
        int next_face_id = wv.E(next_edge_id, 0);
        if (next_face_id == curr_face_id)
        {
            next_face_id = wv.E(next_edge_id, 1);
        }
        if (next_face_id == -1)
        {
            curve.conservativeResize(i+1, 3);
            normal.conservativeResize(i+1, 3);
            break;
        }

        switch (trace_state)
        {
        case GEODESIC:
            curr_dir = mapVectorToAdjacentFace(wv.F, wv.V, wv.edgeVerts,
                next_edge_id, curr_face_id, next_face_id, curr_dir);
            
	    break;
        case FIELD:
	    Eigen::Vector3d parTransport = mapVectorToAdjacentFace(wv.F, wv.V, wv.edgeVerts,
		                    next_edge_id, curr_face_id, next_face_id, curr_dir);
	    Eigen::MatrixXi perm = wv.Ps[next_edge_id];
            if (next_face_id == wv.E(next_edge_id, 1))
            {
                perm.transposeInPlace();
            }
            Eigen::Vector3i curr_vec = Eigen::Vector3i::Zero();
            curr_vec(curr_dir_idx) = 1;
            Eigen::Vector3i next_vec = perm * curr_vec;
            for (int idx = 0; idx < 3; idx++)
            {
                if (abs(next_vec(idx)) > 10e-8)
                {
                    if (next_vec(idx) < 0.)
                    {
                        coeff_dir *= -1.;
                    }
                    curr_dir = wv.Bs[next_face_id] * wv.v(next_face_id, idx);
                    curr_dir_idx = idx;
                    break;
                }
            }
            curr_dir = coeff_dir * curr_dir.normalized() * wv.averageEdgeLength * 1000.;
	    bend(i) = curr_dir.normalized().dot( parTransport.normalized() );
            break;
        }

        curr_face_id = next_face_id;
        curr_edge_id = next_edge_id;
        prev_point = curve.row(i);
        op_v_id = getOpVIdFromEdge(wv, curr_edge_id, curr_face_id);
    }
    curves.push_back(curve);
    normals.push_back(normal);
    modes.push_back(trace_state);
    bending.push_back(bend);
    return;
}

void Trace::save(const std::string &filename)
{
    std::ofstream ofs(filename);
    int ntraces = curves.size();
    ofs << ntraces << std::endl;
    for (int i = 0; i < ntraces; i++)
    {
        int curvesrows = curves[i].rows();
        ofs << curvesrows << std::endl;
        for (int j = 0; j < curvesrows; j++)
        {
            ofs << curves[i](j, 0) << " " << curves[i](j, 1) << " " << curves[i](j, 2) << std::endl;
        }
        int normalsrows = normals[i].rows();
        ofs << normalsrows << std::endl;
        for (int j = 0; j < normalsrows; j++)
        {
            ofs << normals[i](j, 0) << " " << normals[i](j, 1) << " " << normals[i](j, 2) << std::endl;
        }
        int bendingrows = bending[i].rows();
        ofs << bendingrows << std::endl;
        for (int j = 0; j < bendingrows; j++)
        {
            ofs << bending[i][j] << std::endl;
        }

        ofs << (int)modes[i] << std::endl;
    }
}

void Trace::load(const std::string &filename)
{
    std::ifstream ifs(filename);
    if (!ifs)
    {
        std::cerr << "Couldn't load trace file" << std::endl;
        return;
    }
        
    int ntraces;
    ifs >> ntraces;
    curves.resize(ntraces);
    normals.resize(ntraces);
    bending.resize(ntraces);
    modes.resize(ntraces);

    for (int i = 0; i < ntraces; i++)
    {
        int curvesrows;
        ifs >> curvesrows;
        curves[i].resize(curvesrows, 3);
        for (int j = 0; j < curvesrows; j++)
        {
            ifs >> curves[i](j, 0) >> curves[i](j, 1) >> curves[i](j, 2);
        }

        int normalsrows;
        ifs >> normalsrows;
        normals[i].resize(normalsrows, 3);
        for (int j = 0; j < normalsrows; j++)
        {
            ifs >> normals[i](j, 0) >> normals[i](j, 1) >> normals[i](j, 2);
        }

        int bendingrows;
        ifs >> bendingrows;
        bending[i].resize(bendingrows);
        for (int j = 0; j < bendingrows; j++)
        {
            ifs >> bending[i][j];
        }

        int mode;
        ifs >> mode;
        modes[i] = (Trace_Mode)mode;
    }
}
