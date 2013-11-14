/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "Eigen/Core"
#include "MinEig.h"

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www-2.cs.cmu.edu/~ph/869/www/notes/vanishing.txt.
//	
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
    // check
    if (lines.size() < 2)
	{
            fprintf(stderr, "Not enough lines to compute the best fit.");
            abort();
	}

    SVMPoint bestfit;
    list<SVMLine>::const_iterator iter;

    // To accumulate stuff
    typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;

    int numLines = (int) lines.size();
    Matrix3 A = Matrix3::Zero(numLines, 3);	

    // Transformation for numerical stability

    // Note: iterate through the lines list as follows:
    //		for (iter = lines.begin(); iter != lines.end(); iter++) {
    //			...iter is the pointer to the current line...
    //		}
    // Note: Function to find eigenvector with smallest eigenvalue is MinEig(A, eval, evec)
    //
    /******** BEGIN TODO ********/
	int linecount=0;
	for (iter=lines.begin();iter!=lines.end();iter++)
	{
		A(linecount,0)=((iter->pnt1->v)*(iter->pnt2->w))-((iter->pnt1->w)*(iter->pnt2->v));
		A(linecount,1)=((iter->pnt1->w)*(iter->pnt2->u))-((iter->pnt1->u)*(iter->pnt2->w));
		A(linecount,2)=((iter->pnt1->u)*(iter->pnt2->v))-((iter->pnt1->v)*(iter->pnt2->u));
		linecount++;
	}
	Matrix3 M=Matrix3::Zero(3,3);
	for (int i=0;i<numLines;i++)
	{
		M(0,0)+=A(i,0)*A(i,0);M(0,1)+=A(i,0)*A(i,1);M(0,2)+=A(i,0)*A(i,2);
		M(1,0)+=A(i,1)*A(i,0);M(1,1)+=A(i,1)*A(i,1);M(1,2)+=A(i,1)*A(i,2);
		M(2,0)+=A(i,2)*A(i,0);M(2,1)+=A(i,2)*A(i,1);M(2,2)+=A(i,2)*A(i,2);
	}
	double evec[3];
	double eval;
	MinEig(A,eval,evec);
	bestfit.u=evec[0]/evec[2];bestfit.v=evec[1]/evec[2];bestfit.w=evec[2];

    /******** END TODO ********/
	
    return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		plane coordinates. See the following document for more detail.
//		http://www.cs.cornell.edu/courses/cs4670/2012fa/projects/p4/homography.pdf.
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
    int numPoints = points.size();

    /******** BEGIN TODO ********/
	// Choose p,q,r from points to define a unique plane
	Vec4d r = Vec4d(points[0].X, points[0].Y, points[0].Z, points[0].W);
	Vec4d p = Vec4d(points[1].X, points[1].Y, points[1].Z, points[1].W);
	Vec4d q;

	Vec4d rp = p - r;
	rp.normalize();
	
	// Initialize as worst
	double best_angleCOS = 1.0;

	// Loop to find the best q s.t. angle(rp,rq)->90 degrees, i.e. dot(rp/|rp|,rq/|rq|)->0
	for (int cnt_q = 2; cnt_q < numPoints; cnt_q++)
	{
		Vec4d qTmp = Vec4d(points[cnt_q].X, points[cnt_q].Y, points[cnt_q].Z, points[cnt_q].W);
		Vec4d rqTmp = qTmp - r;
		rqTmp.normalize();
		// Compute the angle between rp and rq
		double angleCOS = rp * rqTmp;
		if (abs(angleCOS) < abs(best_angleCOS))
		{
			best_angleCOS = angleCOS;
			q = qTmp;
		}
	}
	
	Vec4d rq = q - r;

	// Normalized already
	Vec4d e_x = rp;
	Vec4d s = (rq * e_x) * e_x;
	Vec4d t = rq - s;
	Vec4d e_y = t;
	e_y.normalize();

	// Ready for u,v
	double u_min = FLT_MAX, v_min = FLT_MAX;
	double u_max = 0, v_max = 0;

	for (int cnt_a = 0; cnt_a < numPoints; cnt_a++)
	{
		Vec4d a = Vec4d(points[cnt_a].X, points[cnt_a].Y, points[cnt_a].Z, points[cnt_a].W);
		Vec4d ra = a - r;
		Vec3d Coord2D = Vec3d(ra*e_x, ra*e_y, 1.0);
		basisPts.push_back(Coord2D);
		if (Coord2D[0] < u_min)
			u_min = Coord2D[0];
		if (Coord2D[0] > u_max)
			u_max = Coord2D[0];
		if (Coord2D[1] < v_min)
			v_min = Coord2D[1];
		if (Coord2D[1] > v_max)
			v_max = Coord2D[1];
	}

	// Output
	uScale = u_max - u_min;
	vScale = v_max - v_min;

    /******** END TODO ********/
}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see
//		http://www.cs.cornell.edu/courses/cs4670/2012fa/projects/p4/homography.pdf.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
    int i;
    int numPoints = (int) points.size();
    assert( numPoints >= 4 );

    basisPts.clear();
    if (isRefPlane) // reference plane
    {
        for (i=0; i < numPoints; i++) {
            Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
            basisPts.push_back(tmp);
        }
    } 
    else // arbitrary polygon
    {
        double uScale, vScale; // unused in this function
        ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
    }

    // A: 2n x 9 matrix where n is the number of points on the plane
    //    as discussed in lecture
    int numRows = 2 * numPoints;
    const int numCols = 9;

    typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
    MatrixType A = MatrixType::Zero(numRows, numCols);

    /******** BEGIN TODO ********/
    /* Fill in the A matrix for the call to MinEig */
printf("TODO: %s:%d\n", __FILE__, __LINE__); 


    double eval, h[9];
    MinEig(A, eval, h);

    H[0][0] = h[0];
    H[0][1] = h[1];
    H[0][2] = h[2];

    H[1][0] = h[3];
    H[1][1] = h[4];
    H[1][2] = h[5];

    H[2][0] = h[6];
    H[2][1] = h[7];
    H[2][2] = h[8];

    /******** END TODO ********/

    // compute inverse of H
    if (H.Determinant() == 0)
        fl_alert("Computed homography matrix is uninvertible \n");
    else
        Hinv = H.Inverse();

    int ii;
    printf("\nH=[\n");
    for (ii=0; ii<3; ii++)
        printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
    printf("]\nHinv=[\n");

    for (ii=0; ii<3; ii++)
        printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);

    printf("]\n\n");
}

