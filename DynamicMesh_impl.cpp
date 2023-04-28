#include "DynamicMesh_impl.h"

int Vertex::nextId=0;

inline bool isLeft(QVector2D a, QVector2D b, QVector2D p){
	return (b.x() - a.x())*(p.y() - a.y()) > (b.y() - a.y())*(p.x() - a.x());
}

bool Vertex::MoveTo(QVector2D newPos, double hardMinHeightFactor, bool checkOnly){
	if(ID()==100){
		Q_ASSERT(true);
	}
	double maxDist=std::numeric_limits<double>::infinity();
	for(int i=0;i<Data().size();i++){
		Triangle * t=(Triangle*)Data()[i];
		for(int j=0;j<3;j++){
			if((*t)(j) == this){ //other two are target edge
				QVector2D one=(*(*t)((j+1)%3));
				QVector2D two=(*(*t)((j+2)%3));
				if(!isLeft(one, two, newPos)){continue;}
				double insetDist=t->OrigMinHeight(j)*hardMinHeightFactor;
				QVector2D insetVec=RotateCW90(two-one).normalized()*insetDist;
				one+=insetVec;
				two+=insetVec;
				QVector2D intersect = ComputeLineIntersection(one, two, (*this), newPos);
				double d=(intersect-(*this)).length();
				//*if we are already on the wrong side of the line (floating point inacc may cause this),
				//dont move unless newPos is on the correct side again or moves it closer to the line
				if(isLeft(one, two, (*this)) && isLeft(one,two, newPos)){
					double newDToline2=(newPos-ProjectPointLine(one,two,newPos)).lengthSquared();
					double dToline2=((*this)-ProjectPointLine(one,two,(*this))).lengthSquared();
					if(newDToline2>dToline2){
						d=0;}
				}//*/
				if(d<maxDist){maxDist=d;}
			}
		}
	}
	double d=((*this) - newPos).length();

	//adjust move if needed
	bool adjusted=false;
	if(d>maxDist){
		qDebug() << "Adjusting illegal move!";
		if(maxDist<1){d=0;}
		else{d=maxDist-1;} //we dont want to put it exactly on the line because floating point inacc may cause probs
		adjusted=true;
	}
	if(checkOnly){return adjusted;}
	if(d!=0){
		//do the actual move
		newPos=(*this)+(newPos-(*this)).normalized()*d;
		setX(newPos.x());
		setY(newPos.y());
	}
	return adjusted;
}
