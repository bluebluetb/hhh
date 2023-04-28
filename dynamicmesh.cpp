#include "dynamicmesh.h"

/*****   TriangleP   *****/
namespace DMesh{
inline QVector2D ProjectPointLine ( QVector2D a , QVector2D b , QVector2D p ) {
	return a + (b - a ) *QVector2D::dotProduct(( p - a ), (b - a )) /QVector2D::dotProduct(( b - a ),(b - a ));
}
}

double TriangleP::MinHeight() const{
	// find smallest height of triangle
	double minHeight = std::numeric_limits<double>::infinity();
	if(HasDuplicates()){return 0;}
	for (int i = 0; i < 3; ++i)
	{
		int j = (i + 1) % 3;
		int k = (j + 1) % 3;
		QVector2D u = (*this)(i);
		QVector2D v = (*this)(j);
		QVector2D w = (*this)(k);
		QVector2D d = ProjectPointLine(u, v, w);
		QVector2D wd = w - d;
		double height = wd.length();
		if (height < minHeight)
		{
			minHeight = height;
		}
	}
	return minHeight;
}

double TriangleP::Area() const{
	double baseLen=(a-b).length();
	double height=(ProjectPointLine(a,b,c)-c).length();
	return baseLen*height/2;
}

bool TriangleP::IsClockwise() const{
	QVector2D centroid = Centroid();

	double dr0 = (*this)(0).x() - centroid.x(),
			dc0 = (*this)(0).y() - centroid.y();
	double dx01 = (*this)(1).x() - (*this)(0).x(),
			dy01 = (*this)(1).y() - (*this)(0).y();

	double df = -dx01 * dc0 + dy01 * dr0;
	return df <= 0;
}

bool MeshVertex::RemoveFromTriangle(MeshTriangle *T){
	for(int i=0;i<Data().size();i++){
		if(Data()[i]==T){
			Data().remove(i);
			return true;
		}
	}
	return false;
}

void MeshTriangle::Set(unsigned int i, MeshVertex* v){
	Q_ASSERT(v!=nullptr);
	MeshVertex * oldV=(*this)(i);
	if(oldV!=nullptr){
		oldV->RemoveFromTriangle(this);
	}
	v->Data().append(this);
	(*this)[i]=v;
}

bool MeshTriangle::Swap(MeshVertex * oldV, MeshVertex * newV){
	for(int i=0;i<3;i++){
		if((*this)(i)==oldV){
			Set(i,newV);
			return true;
		}
	}
	return false;
}

TriangleP MeshTriangle::SwapResult(MeshVertex * oldV, MeshVertex * newV) const{
	TriangleP t(*(*this)(0), *(*this)(1), *(*this)(2));
	for(int i=0;i<3;i++){
		if((*this)(i)==oldV){
			t[i]=*newV;
		}
	}
	return t;
}

void TriangleP::Scale(double factor){
	QVector2D cen=Centroid();
	for(int i=0;i<3;i++){
		(*this)[i]=cen+((*this)[i]-cen)*factor;
	}
}


bool TriangleP::Encloses(const QVector2D &p) const{
	return (SameSide(p,(*this)(0),(*this)(1),(*this)(2)) &&
					SameSide(p,(*this)(1),(*this)(0),(*this)(2)) &&
					SameSide(p,(*this)(2),(*this)(0),(*this)(1)));
}
