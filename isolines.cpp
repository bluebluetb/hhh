#include "isolines.h"

QPair<double,double> LinearInterp(QPair<double, QPair<double,double> > p1,
																	QPair<double, QPair<double,double> > p2,
																	float value)
{
	QPair<double,double> p;
	double P=p2.second.first;
	double Q=p1.second.first;//TODO: netjes maken
	if(p1.first != p2.first){
		p.first = Q + (P - Q)*((value - p1.first)/(p2.first - p1.first));
		p.second = p1.second.second + (p2.second.second - p1.second.second)*((value - p1.first)/(p2.first - p1.first));
	}
	else{
		p = qMakePair(p1.second.first,p1.second.second);
	}
	return p;
}

QPair<double,double> LinearInterp(QPair<double, QPair<int,int> > p1,
																	QPair<double, QPair<int,int> > p2,
																	float value){
	return LinearInterp(qMakePair(p1.first,qMakePair(double(p1.second.first),double(p1.second.second))),
											qMakePair(p2.first,qMakePair(double(p2.second.first),double(p2.second.second))),
											value
											);
}

const short int edgeTable[16]={
	0b00000000, 0b00001001, 0b00000011, 0b00001010,
	0b00000110, 0b10001111, 0b00000101, 0b00001100,
	0b00001100, 0b00000101, 0b10001111, 0b00000110,
	0b00001010, 0b00000011, 0b00001001, 0b00000000
};

const short int lineTable[32][10]={
	{-1,-1,-1,-1,-1,-1},
	{ 0, 3,-1,-1,-1,-1},
	{ 0, 1,-1,-1,-1,-1},
	{ 1, 3,-1,-1,-1,-1},
	{ 1, 2,-1,-1,-1,-1},
	{ 0, 1, 2, 3,-1,-1},
	{ 0, 2,-1,-1,-1,-1},
	{ 2, 3,-1,-1,-1,-1},
	{ 2, 3,-1,-1,-1,-1},
	{ 0, 2,-1,-1,-1,-1},
	{ 0, 3, 2, 1,-1,-1},
	{ 1, 2,-1,-1,-1,-1},
	{ 1, 3,-1,-1,-1,-1},
	{ 0, 1,-1,-1,-1,-1},
	{ 0, 3,-1,-1,-1,-1},
	{-1,-1,-1,-1,-1,-1},

	{-1,-1,-1,-1,-1,-1},
	{ 0, 3,-1,-1,-1,-1},
	{ 0, 1,-1,-1,-1,-1},
	{ 1, 3,-1,-1,-1,-1},
	{ 1, 2,-1,-1,-1,-1},
	{ 0, 3, 2, 1,-1,-1},
	{ 0, 2,-1,-1,-1,-1},
	{ 2, 3,-1,-1,-1,-1},
	{ 2, 3,-1,-1,-1,-1},
	{ 0, 2,-1,-1,-1,-1},
	{ 0, 1, 2, 3,-1,-1},
	{ 1, 2,-1,-1,-1,-1},
	{ 1, 3,-1,-1,-1,-1},
	{ 0, 1,-1,-1,-1,-1},
	{ 0, 3,-1,-1,-1,-1},
	{-1,-1,-1,-1,-1,-1}
};

Isolines::Isolines(QRectF window, int gridSize, int nLines): window(window), gridSize(gridSize), nLines(nLines)
{
}

void Isolines::DrawIsoline(float value){
	glBegin(GL_LINES);
	//glColor3f(1,0,0);
	double wn = window.width() / (double) (gridSize + 1);   // Grid cell width
	double hn = window.height() / (double) (gridSize + 1);  // Grid cell heigh
	for (int y = 0; y < gridSize - 1; y++)
	{
		for (int x = 0; x < gridSize - 1; x++) {;
			//vertex initialization
			QPair<double, QPair<int, int> > verts[4]; //4 corners of the current cell
			verts[0].first = GetScalarValue(QVector2D(x*wn+window.left(),y*hn+window.top()));
			verts[0].second = qMakePair(x, y);
			verts[1].first = GetScalarValue(QVector2D((x+1)*wn+window.left(),y*hn+window.top()));
			verts[1].second = qMakePair(x+1, y);
			verts[2].first = GetScalarValue(QVector2D((x+1)*wn+window.left(),(y+1)*hn+window.top()));
			verts[2].second = qMakePair(x+1, y+1);
			verts[3].first = GetScalarValue(QVector2D(x*wn+window.left(),(y+1)*hn+window.top()));
			verts[3].second = qMakePair(x, y+1);

			//which corners are below/equal to our value
			short squareIndex = 0;
			for(int n=0; n < 4; n++){
				if(verts[n].first <= value) squareIndex |= 1 << (n);
			}
			//get interpolated values
			QPair<double,double>intVerts[4];
			if(!edgeTable[squareIndex]) continue;
			if(edgeTable[squareIndex] & 1) intVerts[0] = LinearInterp(verts[0], verts[1], value);
			if(edgeTable[squareIndex] & 2) intVerts[1] = LinearInterp(verts[1], verts[2], value);
			if(edgeTable[squareIndex] & 4) intVerts[2] = LinearInterp(verts[2], verts[3], value);
			if(edgeTable[squareIndex] & 8) intVerts[3] = LinearInterp(verts[3], verts[0], value);
			if(edgeTable[squareIndex] & 128){
				//calculate assymtotic deZider
				float D=verts[0].first;
				float B=verts[1].first-verts[0].first;
				float C=verts[3].first-verts[0].first;
				float A=verts[0].first+verts[2].first-verts[1].first-verts[3].first;
				if(A==0 || (D-(B*C/A))>=value){
					squareIndex|=0b00010000;
				}
			}
			for (int n=0; lineTable[squareIndex][n] != -1; n+=2) {
				float f=wn + intVerts[lineTable[squareIndex][n+1]].first * wn;
				float s=hn + intVerts[lineTable[squareIndex][n+1]].second * hn;
				float f2=wn + intVerts[lineTable[squareIndex][n]].first * wn;
				float s2=hn + intVerts[lineTable[squareIndex][n]].second * hn;

				glVertex2d(f+window.left(),s+window.top());
				glVertex2d(f2+window.left(),s2+window.top());
			}
		}
	}
	glEnd();

}


void Isolines::DrawIsolines() {
	//return;
	float minimum=INT_MAX;
	float maximum=INT_MIN;

	for(int i=0;i<gridSize;i++){
		for(int j=0;j<gridSize;j++){
			double x=(window.width()/double(gridSize-1))*i+window.left();
			double y=(window.height()/double(gridSize-1))*i+window.top();

			minimum=min((double)minimum,GetScalarValue(QVector2D(x,y)));
			maximum=max((double)maximum,GetScalarValue(QVector2D(x,y)));
		}
	}
	float val;
	//DrawIsoline(2);
	//DrawIsoline(7);

	//return;
	for(float i=1;i<=nLines;i+=1){
		val=minimum+((maximum-minimum)/(nLines+1))*i;
		DrawIsoline(val);
		qDebug()<<"drawing isoline "<<val;
	}
}

void DisIsolines::UpdateVals(){
	realPoints.clear();
	VertexIterator it=p->mesh->VBegin();
	for (; it!=p->mesh->VEnd(); it.Next()) {
		QVector2D v=it.GetV();
		if(v.x()<window.left()){window.setLeft(v.x());}
		if(v.x()>window.right()){window.setRight(v.x());}
		if(v.y()<window.bottom()){window.setBottom(v.y());}
		if(v.y()>window.top()){window.setTop(v.y());}
	}
	//retriangulate real datapoints for interpolation later on
	triangulation.clear();
	std::vector< Delaunay::Point > v;
	for (it=p->mesh->VBegin(); it!=p->mesh->VEnd(); it.Next()) {
		if(it.Vertex().IsReal()){
			Delaunay::Point tempP;
			tempP[0] = it.x();
			tempP[1] = it.y();
			v.push_back(tempP);
			realPoints.append(it.VertexPointer());
		}
	}
	Delaunay delobject(v);
	delobject.Triangulate();
	for (Delaunay::fIterator fit = delobject.fbegin(); fit != delobject.fend(); ++fit) {
		triangulation.append(DMesh::Triangle(realPoints[delobject.Org(fit)], realPoints[delobject.Dest(fit)], realPoints[delobject.Apex(fit)]));
	}
}

void DisIsolines::DrawIsolines(){
	UpdateVals();
	Isolines::DrawIsolines();
}

//QVector2D ProjectPointLine ( QVector2D a , QVector2D b , QVector2D c ) {
//	return a + (b - a ) *QVector2D::dotProduct(( c - a ), (b - a )) /QVector2D::dotProduct(( b - a ),(b - a ));
//}

QVector2D ProjectPointSegment ( QVector2D a , QVector2D b , QVector2D c ) {
	double r = QVector2D::dotProduct((b - a ),(b - a ));
	if ( fabs ( r ) < 1e-8 ) return a ;
	r = QVector2D::dotProduct(( c - a ),(b - a )) / r ;
	if ( r < 0) return a ;
	if ( r > 1) return b ;
	return a + (b - a ) * r ;
}



double DisIsolines::GetScalarValue(QVector2D pos){ //assumes atleast 2 real points
	//at least 2 real points
	/*
	double minDist=9999999;
	int min=0;
	for(int i=offset;i<p->size();i++){
		double d=(pos-p->X.GetV(i)).length();
		if(d<minDist){
			minDist=d;
			min=i;
		}
	}
	return p->X.VertexData(min).value;
	*/
	//find triangle in which pos resides
	int triangle=-1;
	for(int i=0;i<triangulation.size();i++){
		if(triangulation[i].toTriangleP().Contains(pos)){triangle=i;break;}
	}
	if(triangle>=0){
		//point inside a triangle
		QVector3D P=*(triangulation[triangle](0));
		QVector3D Q=*(triangulation[triangle](1));
		QVector3D R=*(triangulation[triangle](2));
		P.setZ(((Vertex*)triangulation[triangle](0))->value);
		Q.setZ(((Vertex*)triangulation[triangle](1))->value);
		R.setZ(((Vertex*)triangulation[triangle](2))->value);


		double a=-(R.y()*Q.z()-P.y()*Q.z()-R.y()*P.z()+P.z()*Q.y()+R.z()*P.y()-Q.y()*R.z());
		double b= (P.y()*R.x()+Q.y()*P.x()+R.y()*Q.x()-Q.y()*R.x()-P.y()*Q.x()-R.y()*P.x());
		double c= (Q.z()*R.x()+P.z()*Q.x()+R.z()*P.x()-P.z()*R.x()-Q.z()*P.x()-Q.x()*R.z());
		double d=-a*P.x()-b*P.z()-c*P.y();

		return -(a*pos.x()+c*pos.y()+d)/b; //z of point
	}
	else{
		//point outside of all triangles, project onto nearest edge and linear interpolate
		double minDist=std::numeric_limits<double>::infinity(), minDist2=std::numeric_limits<double>::infinity();
		int min=-1, min2=-1;
		for(int i=0;i<realPoints.size();i++){
			double d=(pos-(*realPoints[i])).length();
			if(d<minDist){minDist2=minDist;minDist=d;min2=min;min=i;}
			else{
				if(d<minDist2){minDist2=d;min2=i;}
			}
		}
		if(min==-1 || min2==-1){return 0;}
		QPair<QVector2D, QVector2D> segment=qMakePair(QVector2D(*realPoints[min]),QVector2D(*realPoints[min2]));
		QVector2D projected=ProjectPointSegment(segment.first,segment.second,pos);
		if(projected==segment.first){return realPoints[min]->value;}
		if(projected==segment.second){return realPoints[min2]->value;}
		double segLength=(segment.first-segment.second).length();
		double ratio=(segment.first-projected).length()/segLength;
		return ratio*realPoints[min2]->value+(1-ratio)*realPoints[min]->value;
	}
}
