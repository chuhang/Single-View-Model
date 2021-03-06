/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
void ImgView::sameXY()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	if( refPointOffPlane == NULL )
	{
		fl_alert("Need to specify the reference height first.");
		return;
	}

	/******** BEGIN TODO ********/

	// See the lecture note on measuring heights
	// using a known point directly below the new point.
	//Mat3d matHinv=Mat3d(Hinv[0][0],Hinv[0][1],Hinv[0][2],Hinv[1][0],Hinv[1][1],Hinv[1][2],Hinv[2][0],Hinv[2][1],Hinv[2][2]);
    //Mat3d matH=Mat3d(H[0][0],H[0][1],H[0][2],H[1][0],H[1][1],H[1][2],H[2][0],H[2][1],H[2][2]);

	Vec3d t0=Vec3d(newPoint.u/newPoint.w,newPoint.v/newPoint.w,1);
	Vec3d r=Vec3d(refPointOffPlane->u/refPointOffPlane->w,refPointOffPlane->v/refPointOffPlane->w,1);
	Vec3d vx=Vec3d(xVanish.u,xVanish.v,xVanish.w);
    Vec3d vy=Vec3d(yVanish.u,yVanish.v,yVanish.w);
    Vec3d vz=Vec3d(zVanish.u,zVanish.v,zVanish.w);
	Mat3d matH=Mat3d(H[0][0],H[0][1],H[0][2],H[1][0],H[1][1],H[1][2],H[2][0],H[2][1],H[2][2]);

	Vec3d horizon=cross(vx,vy);
	Vec3d b0=matH*Vec3d(knownPoint.X,knownPoint.Y,1);
	Vec3d b=matH*Vec3d(refPointOffPlane->X,refPointOffPlane->Y,1);
	Vec3d v=cross(cross(b,b0),horizon);
	printf( "the v: (%e, %e, %e)\n", v[0], v[1], v[2] );
	Vec3d t;
	if (v[0]==0 && v[1]==0 && v[2]==0)
	{
		t=t0-b0+b;
	}
	else
	{
		Vec3d rb=cross(r,b);
		t=cross(cross(v,t0),rb);
	}
	printf( "the t: (%e, %e, %e)\n", t[0], b[1], r[2] );

	//printf( "the weights: (%e, %e, %e)\n", t[2], b[2], r[2] );
	t[0]=t[0]/t[2];t[1]=t[1]/t[2];t[2]=1;
	b[0]=b[0]/b[2];b[1]=b[1]/b[2];b[2]=1;
	vz[0]=vz[0]/vz[2];vz[1]=vz[1]/vz[2];vz[2]=1;
	r[0]=r[0]/r[2];r[1]=r[1]/r[2];r[2]=1;
	
	double crossratio=sqrt((t-b)*(t-b)) * sqrt((vz-r)*(vz-r)) / sqrt((r-b)*(r-b)) / sqrt((vz-t)*(vz-t));
	double h=crossratio*referenceHeight;
	Vec3d b0t0 = t0 - b0;
	Vec3d bvz = vz - b;
	double dotprdt = b0t0[0]*bvz[0] + b0t0[1]*bvz[1] + b0t0[2]*bvz[2];
	if (dotprdt < 0) h *= -1;

	newPoint.X=knownPoint.X;
	newPoint.Y=knownPoint.Y;
	newPoint.Z=h;
	newPoint.W=1;
	printf("sameXY!\n");
	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	/******** BEGIN TODO ********/
	Mat3d matHinv=Mat3d(Hinv[0][0],Hinv[0][1],Hinv[0][2],Hinv[1][0],Hinv[1][1],Hinv[1][2],Hinv[2][0],Hinv[2][1],Hinv[2][2]);
    Mat3d matH=Mat3d(H[0][0],H[0][1],H[0][2],H[1][0],H[1][1],H[1][2],H[2][0],H[2][1],H[2][2]);

    Vec3d t1=Vec3d(newPoint.u,newPoint.v,newPoint.w);
    Vec3d m0=Vec3d(knownPoint.u,knownPoint.v,knownPoint.w);
    Vec3d b1=Vec3d(newPoint.u,newPoint.v,newPoint.w);
	Vec3d nP=Vec3d(newPoint.u,newPoint.v,newPoint.w);

    if (knownPoint.Z != 0)
	{
        Vec3d vx=Vec3d(xVanish.u,xVanish.v,xVanish.w);
        Vec3d vy=Vec3d(yVanish.u,yVanish.v,yVanish.w);
        Vec3d vz=Vec3d(zVanish.u,zVanish.v,zVanish.w);      
        // vertical line passing through the newpoint
		Vec3d t1b1=cross(t1,vz);
		// horizon
        Vec3d horizon=cross(vx,vy);
		// vanishing point
        Vec3d v=cross(cross(t1,m0),horizon);               
        Vec3d b0=matH*Vec3d(knownPoint.X,knownPoint.Y,1);
        // if t1m0 parallel to horizon
		if (v[0]==0 && v[1]==0 && v[2]==0)
			b1 = cross(b0,t1b1);
		else
			b1=cross(cross(b0,v),t1b1);
    }

    b1[0]=b1[0]/b1[2];
    b1[1]=b1[1]/b1[2];
    b1[2]=b1[2]/b1[2];    
    nP=matHinv*b1;     
    newPoint.X=nP[0]/nP[2];
    newPoint.Y=nP[1]/nP[2];
    newPoint.Z=knownPoint.Z;
    newPoint.W=1;
	printf("sameZPlane!\n");
	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}


