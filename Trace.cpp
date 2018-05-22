#if 0

#include "Trace.h"
#include "VectorUtils.h"
#include "Distance.h"
#include "Weave.h"
#include "Surface.h"

#include <fstream>
#include <iostream>

#include <igl/point_mesh_squared_distance.h>

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
    
    int counter = 0;
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


     
        if (npoints > 80 && npoints < 250 && counter < 50 || true)
        {
            counter++;
            curves.push_back(curve);
            normals.push_back(normal);
            bending.push_back(bend);
            modes.push_back(Trace_Mode::FIELD); // Make another render option
        }
    }
    curve_file.close();
    normal_file.close();
}


void Trace::loadGeneratedCurves(std::vector<std::vector<Eigen::Vector3d> > isoLines, 
                                std::vector<std::vector<Eigen::Vector3d> > isoNormal)
{

    curves.clear();
    normals.clear();
    bending.clear();
    modes.clear();
    for (int i = 0; i < isoLines.size(); i++)
    {
        double npoints = isoLines[i].size();
        Eigen::MatrixXd curve = Eigen::MatrixXd::Zero(npoints, 3);
        Eigen::MatrixXd normal = Eigen::MatrixXd::Zero(npoints, 3);
        Eigen::VectorXd bend = Eigen::VectorXd::Zero(npoints); // Not implemented yet
        if ( npoints > 3 )
        {
            for (int j = 0; j < npoints; j++)
            {
                curve.row(j) = isoLines[i][j];
                normal.row(j) = isoNormal[i][j];

            }
            curves.push_back(curve);
            normals.push_back(normal);
            bending.push_back(bend);
            modes.push_back(Trace_Mode::FIELD);
        }
    }
}


double angle(Eigen::Vector3d v, Eigen::Vector3d w, Eigen::Vector3d n)
{
    return 2.0 * atan2( v.cross(w).dot(n), v.norm() * w.norm() + v.dot(w));
}


void Trace::logRibbonsToFile(std::string foldername, std::string filename, const Weave &wv)
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




    // Project along a geodesic curve at the end point of each rod that ends at a singularity.  
    
    Eigen::VectorXd sqrD;
    Eigen::VectorXi nearFace;
    Eigen::MatrixXd C;
    Eigen::MatrixXd point = Eigen::MatrixXd::Zero(1,3); 
    igl::AABB<Eigen::MatrixXd,3> tree; 
   
    tree.init(wv.fs->data().V,wv.fs->data().F);

 
 // need to deal with boundary case
    int stepstoextend = 0;
    int backoff = 0; // The very ends of curves seem generally bad.  Taking a few steps to back off

    for (int i = 0; i < curves.size(); i++)
    { 
        int len = curves[i].rows();
        int cols = curves[i].cols();
        curves[i].conservativeResize(len + stepstoextend - backoff, cols);
        normals[i].conservativeResize(len + stepstoextend - backoff, cols);

        Eigen::Vector3d endpoint = curves[i].row(len-1 - backoff);
        Eigen::Vector3d prevpoint = curves[i].row(len-2 - backoff);
        Eigen::Vector3d curr_dir = endpoint - prevpoint;
        // hack-y way of figuring out the edge of the end-point 
        point.row(0) = prevpoint + curr_dir * .01;

        curr_dir *= 1000.0; // Assumes mesh is roughly deluanay

        tree.squared_distance(wv.fs->data().V,wv.fs->data().F,point,sqrD,nearFace,C);
        int curr_face_id = nearFace(0);

        TracePoint tp;
        startTraceFromPoint(wv, curr_face_id, point.row(0), curr_dir, tp);
//        std::cout << point.row(0) - prevpoint << " should be same \n";
        int curr_edge_id = tp.edge_id;
        getNextTracePoint(wv, curr_face_id, tp.edge_id, prevpoint, tp.op_v_id, curr_dir, tp); 
//        std::cout << tp.point - endpoint << " endpoints should be same \n";

        for (int j = 0; j < stepstoextend; j++)
        {
            curr_dir = mapVectorToAdjacentFace(wv.fs->data().F, wv.fs->data().V, wv.fs->data().edgeVerts,
                                               tp.edge_id, curr_face_id, tp.face_id, curr_dir);
            curr_face_id = tp.face_id;  
            curves[i].row(len + j - backoff)  = tp.point;
            normals[i].row(len + j - backoff) = tp.n;
            std::cout << normals[i] << "\n";

            getNextTracePoint(wv, curr_face_id, tp.edge_id, tp.point, tp.op_v_id, curr_dir, tp); 
        }

        Eigen::MatrixXd revcurve = curves[i];
        Eigen::MatrixXd revnorm = normals[i];
        int c_len = curves[i].rows();
        for (int r = 0; r < c_len; r++)
        {
            revcurve.row(r) = curves[i].row(c_len - r - 1);
            revnorm.row(r) = normals[i].row(c_len - r - 1);
        }
        curves[i] = revcurve;
        normals[i] = revnorm;
    }   

    // Cut strips at areas of high curvature.

    std::vector<Eigen::MatrixXd> splitcurves;
    std::vector<Eigen::MatrixXd> splitnormals;

    double maxcurvature = 2000; // angle - currently disabled
    double minrodlen = 2;

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









    // Find collisions between rods
    std::vector<Collision> collisions;
    for (int i = 0; i < splitcurves.size(); i++)
    {
        for (int j = i; j < splitcurves.size(); j++)
        {
            computeIntersections(i, j, splitcurves, collisions);
        }
    }

    // Decimate and log rods
    std::vector<Eigen::VectorXd> desc_maps;
    std::vector< std::vector<Eigen::Vector3d> > desc_curves; //eeew
    std::vector< std::vector<Eigen::Vector3d> > desc_normals;

    double decimation_factor = 1.; // make this bigger and add intelligent subdivision

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
        //        max_length = 0.;
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
            if (seg_length > max_length / decimation_factor)
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
            myfile << " 0.02";
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
        if(wv.fs->data().faceNeighbors(faceId, j) == -1)
            continue;
	Eigen::VectorXi e = wv.fs->data().edgeVerts.row(wv.fs->data().faceNeighbors(faceId, j));
	Eigen::Vector3d e_test =  wv.fs->data().V.row( e(0) ) - wv.fs->data().V.row( e(1) );
	e_test.normalize();
        Eigen::Vector3d point_test = wv.fs->data().V.row( e(0) );
        point_test = ( prev_point - point_test );
        point_test.normalize();
	if ( e_test.dot( point_test ) > .99 )
	{
	    Eigen::VectorXi e_next = wv.fs->data().edgeVerts.row(wv.fs->data().faceNeighbors(faceId, (j + 1) % 3 ));
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
    Eigen::Vector2i e = wv.fs->data().edgeVerts.row(curr_edge);
    for (int i = 0; i < 3; i++)
    {
        if (wv.fs->data().F(faceId, i) != e(0) && wv.fs->data().F(faceId, i) != e(1))
        {
            return wv.fs->data().F(faceId, i);
        }
    }
    assert(false);
    return -1;
}





void Trace::getNextTracePoint(const Weave &wv, 
                              int curr_face_id, int curr_edge_id, 
                              const Eigen::Vector3d prev_point, int op_v_id,  
                              const Eigen::Vector3d curr_dir, TracePoint &nextTrace)
{

        
}


/*
void Trace::traceKSteps(const Weave &wv, int traceIdx, int faceId, int steps)
{


}
*/


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

#endif