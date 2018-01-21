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

    double npoints_d, dum1, dum2;
    double dum3, dum4, dum5;
    while (normal_file >> dum3 >> dum4 >> dum5)
    {
        curve_file >> npoints_d >> dum1 >> dum2;
        double npoints = npoints_d;
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

        std::cout << normal.row(npoints - 1) << " " << npoints;



        if (npoints > 20)
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

double angle(Eigen::Vector3d v, Eigen::Vector3d w, Eigen::Vector3d n)
{
    return 2.0 * atan2( v.cross(w).dot(n), v.norm() * w.norm() + v.dot(w));
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


    // Cut strips at areas of high curvature.

    std::vector<Eigen::MatrixXd> splitcurves;
    std::vector<Eigen::MatrixXd> splitnormals;

    double maxcurvature = 30;
    double minrodlen = 30;

    for(int i=0; i<curves.size(); i++)
    {
        if(curves[i].rows() < 3)
            continue;
        std::vector<Eigen::Vector3d> curpts;
        std::vector<Eigen::Vector3d> curnormals;
        curpts.push_back(curves[i].row(0));
        curpts.push_back(curves[i].row(1));
        curnormals.push_back(normals[i].row(0));
        curnormals.push_back(normals[i].row(1));
        Eigen::Vector3d prevedge = curves[i].row(1) - curves[i].row(0);
        for(int j=2; j<curves[i].rows(); j++)
        {
            Eigen::Vector3d nextedge = curves[i].row(j).transpose() - curpts.back();
            Eigen::Vector3d prevproj = prevedge - prevedge.dot(normals[i].row(j)) * normals[i].row(j).transpose();
            if(fabs(angle(prevproj, nextedge, normals[i].row(j)))/(prevedge.norm() + nextedge.norm()) > maxcurvature)
            {
                // cut
                if(curpts.size() > minrodlen)
                {
                    Eigen::MatrixXd newcurve(curpts.size(), 3);
                    Eigen::MatrixXd newnormal(curpts.size(), 3);
                    for(int j=0; j<curpts.size(); j++)
                    {
                        newcurve.row(j) = curpts[j].transpose();
                        newnormal.row(j) = curnormals[j].transpose();
                    }
                    splitcurves.push_back(newcurve);
                    splitnormals.push_back(newnormal);
                }
                curpts.clear();
                curnormals.clear();
                curpts.push_back(curves[i].row(j-1));
                curpts.push_back(curves[i].row(j));
                curnormals.push_back(normals[i].row(j-1));
                curnormals.push_back(normals[i].row(j));
            }
            else
            {
                curpts.push_back(curves[i].row(j));
                curnormals.push_back(normals[i].row(j));
            }
            prevedge = nextedge;
        }
        if(curpts.size() > minrodlen)
        {
            Eigen::MatrixXd newcurve(curpts.size(), 3);
            Eigen::MatrixXd newnormal(curpts.size(), 3);
            for(int j=0; j<curpts.size(); j++)
            {
                newcurve.row(j) = curpts[j].transpose();
                newnormal.row(j) = curnormals[j].transpose();
            }
            splitcurves.push_back(newcurve);
            splitnormals.push_back(newnormal);
        }             
    }


    // Project along a geodesic curve at the end point of each rod that ends at a singularity.  



    // Split each segment into k peices




    // Find collisions between rods
    std::vector<Collision> collisions;
    for (int i = 0; i < splitcurves.size(); i++)
    {
        for (int j = i; j < splitcurves.size(); j++)
        {
            computeIntersections(i, j, collisions, splitcurves);
        }
    }

    // Decimate and log rods
    std::vector<Eigen::VectorXd> desc_maps;
    std::vector< std::vector<Eigen::Vector3d> > desc_curves; //eeew
    std::vector< std::vector<Eigen::Vector3d> > desc_normals;
    for (int curveId = 0; curveId < splitcurves.size(); curveId++)
    {
        Eigen::MatrixXd curve = splitcurves[curveId];
        Eigen::MatrixXd curveNormals = splitnormals[curveId];

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
                prev_point = cnew.back();
            }
        }

        desc_maps.push_back(desc_mapping);
        desc_curves.push_back(cnew);
        desc_normals.push_back(nnew);
    }


    std::vector<Collision> desc_collisions;
    for (int i = 0; i < collisions.size(); i++)
    {
        Collision col = collisions[i];
        std::vector<Eigen::Vector3d> &c1 = desc_curves[col.rod1];
        std::vector<Eigen::Vector3d> &c2 = desc_curves[col.rod2];
        int idx1 = desc_maps[col.rod1](col.seg1);
        int idx2 = desc_maps[col.rod2](col.seg2);
        if(idx1 < c1.size() - 1 && idx2 < c2.size() - 1)
        {
            Collision newcol(col.rod1, col.rod2, idx1, idx2);
            desc_collisions.push_back(newcol);
        }
    }

    // Write Header 
    myfile << desc_curves.size() << std::endl;;
    myfile << desc_collisions.size() << std::endl;;
  //  myfile << 0 << std::endl;;
    myfile << "0.001"  << std::endl;;
    myfile << "1e+08"  << std::endl;;
    myfile << "1"  << std::endl  << std::endl  << std::endl;;


    for (int curveId = 0; curveId < desc_curves.size(); curveId++)
    {
        std::vector<Eigen::Vector3d> &cnew = desc_curves[curveId];
        std::vector<Eigen::Vector3d> &nnew = desc_normals[curveId];
        myfile << cnew.size() << "\n";
        myfile << "0\n";
        for (int i = 0; i < cnew.size(); i++)
        {
            myfile << cnew[i](0) << " " << cnew[i](1) << " " << cnew[i](2) << " ";
        }
        myfile << "\n";

        for (int i = 0; i < nnew.size(); i++)
        {
            myfile << nnew[i](0) << " " << nnew[i](1) << " " << nnew[i](2) << " ";
        }
        myfile << "\n";

        for (int i = 0; i < cnew.size()-1; i++)
        {
            myfile << " 0.005";
        }
        myfile << "\n";

    }
    for (int i = 0; i < desc_collisions.size(); i++)
    {
        Collision col = desc_collisions[i];
        std::vector<Eigen::Vector3d> c1 = desc_curves[col.rod1];
        std::vector<Eigen::Vector3d> c2 = desc_curves[col.rod2];
        int idx1 = col.seg1;
        int idx2 = col.seg2;

        double p0bary, p1bary, q0bary, q1bary;

        if(idx1 >= c1.size()-1) exit(-1);

        Eigen::Vector3d dist = Distance::edgeEdgeDistance(c1[idx1],
            c1[idx1 + 1],
            c2[idx2],
            c2[idx2 + 1],
            p0bary, p1bary, q0bary, q1bary);

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

void Trace::computeIntersections(int curveIdx1, int curveIdx2, std::vector<Collision> &collisions, std::vector<Eigen::MatrixXd> &splitcurves)
{
    Eigen::MatrixXd c1 = splitcurves[curveIdx1];
    Eigen::MatrixXd c2 = splitcurves[curveIdx2];

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


void Trace::startTraceFromPoint(const Weave &wv,
                                int curr_face_id, 
                                const Eigen::Vector3d startpoint, 
                                const Eigen::Vector3d curr_dir, TracePoint &startPoint)
{

    // Project backwards to initialize.  Kinda hacky.
    double min_dist = std::numeric_limits<double>::infinity();

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

    startPoint.edge_id = curr_edge_id;
    startPoint.point = prev_point;

}

void Trace::getNextTracePoint(const Weave &wv, 
                              int curr_face_id, int curr_edge_id, 
                              const Eigen::Vector3d prev_point, int op_v_id,  
                              const Eigen::Vector3d curr_dir, TracePoint &nextTrace)
{

        Eigen::Vector3d op_vertex = wv.V.row(op_v_id);
        Eigen::Vector3d split = op_vertex - prev_point;

        double split_len = split.norm();
        split = split / split_len;
        Eigen::Vector3d n = faceNormal(wv.F, wv.V, curr_face_id);
        Eigen::Vector3d perp = split.cross(n);
        nextTrace.n = n;

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
            nextTrace.edge_id = -1;
            nextTrace.face_id = -1;
            return;
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

        int next_face_id = wv.E(next_edge_id, 0);
        if (next_face_id == curr_face_id)
        {
            next_face_id = wv.E(next_edge_id, 1);
        }
        nextTrace.face_id = next_face_id;
        nextTrace.edge_id = next_edge_id;
        nextTrace.point = (op_v1 * q0bary + op_v2 * q1bary);
}


/*
void Trace::traceKSteps(const Weave &wv, int traceIdx, int faceId, int steps)
{


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

    Eigen::Vector3d startpoint;
    startpoint.setZero();
    for (int i = 0; i < 3; i++)
    {
        startpoint += wv.V.row(wv.F(faceId, i));
    }
    startpoint /= 3.0;


    TracePoint tp;
    startTraceFromPoint(wv, curr_face_id, startpoint, curr_dir, tp); 
    int curr_edge_id = tp.edge_id;
    Eigen::Vector3d prev_point = tp.point;

    assert(curr_edge_id > -1);
    int op_v_id = getOpVIdFromEdge(wv, curr_edge_id, curr_face_id);
    Eigen::VectorXd bend = Eigen::VectorXd::Zero(steps);
    for (int i = 1; i < steps; i++)
    {
            
        getNextTracePoint(wv, curr_face_id, curr_edge_id, prev_point, op_v_id, curr_dir, tp);  
        curve.row(i) = tp.point;
        normal.row(i) = tp.n;


        // stop if we hit vertex
        if (tp.edge_id == -1)
        {
            curve.conservativeResize(i, 3);
            normal.conservativeResize(i, 3);
            break;
        }
        // stop if we hit boundary
        if (tp.face_id == -1)
        {
            curve.conservativeResize(i+1, 3);
            normal.conservativeResize(i+1, 3);
            break;
        }

        switch (trace_state)
        {
        case GEODESIC:
            curr_dir = mapVectorToAdjacentFace(wv.F, wv.V, wv.edgeVerts,
                           tp.edge_id, curr_face_id, tp.face_id, curr_dir);
            
	    break;
        case FIELD:
	        Eigen::Vector3d parTransport = mapVectorToAdjacentFace(wv.F, wv.V, wv.edgeVerts,
		                                   tp.edge_id, curr_face_id, tp.face_id, curr_dir);
	        Eigen::MatrixXi perm = wv.Ps[tp.edge_id];
            if (tp.face_id == wv.E(tp.edge_id, 1))
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
                    curr_dir = wv.Bs[tp.face_id] * wv.v(tp.face_id, idx);
                    curr_dir_idx = idx;
                    break;
                }
            }
            curr_dir = coeff_dir * curr_dir.normalized() * wv.averageEdgeLength * 1000.;
	        bend(i) = curr_dir.normalized().dot( parTransport.normalized() );
            break;
        }

        curr_face_id = tp.face_id;
        curr_edge_id = tp.edge_id;
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
