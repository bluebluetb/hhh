#include "geom.h"

inline double Cross(QVector2D a, QVector2D b){
	return a.x()*b.y()-a.y()*b.x();
}

bool LinesParallel ( QVector2D a , QVector2D b , QVector2D c , QVector2D d ) {
	return fabs (Cross(b-a,c-d)) < EPS ;
}


bool LinesCollinear ( QVector2D a , QVector2D b
											, QVector2D c , QVector2D d ) {
	return LinesParallel (a , b , c , d )
			&& fabs (Cross(a-b,a-c)) <
			EPS
			&& fabs (Cross(c-d,c-a)) <
			EPS ;
}

inline double dist2(QVector2D a, QVector2D b){
	return QVector2D::dotProduct(a,b);
}

bool SegmentsIntersect( QVector2D a , QVector2D b , QVector2D c , QVector2D d ) {
	if ( LinesCollinear (a , b , c , d ) ) {
		if ( dist2 (a , c ) < EPS || dist2 (a , d ) < EPS ||
				 dist2 (b , c ) < EPS || dist2 (b , d ) < EPS ) return true ;
		if (QVector2D::dotProduct(c-a,c-b) > 0 && QVector2D::dotProduct(d-a,d-b) > 0 && QVector2D::dotProduct(c-b,d-b) > 0)
			return false ;
		return true ;
	}
	if (Cross(d-a, b-a) * Cross(c-a, b-a) > 0) return false ;
	if (Cross(a-c, d-c) * Cross(b-c, d-c) > 0) return false ;
	return true ;
}

QVector2D ComputeLineIntersection( QVector2D a , QVector2D b , QVector2D c , QVector2D d ) {
	b=b-a;
	d=c-d;
	c=c-a;
	return a+b*Cross(c,d)/Cross(b,d) ;
}


bool SegmentLineIntersect(QVector2D p1, QVector2D p2, QVector2D a, QVector2D b){
	QVector3D cp1 = QVector3D::crossProduct((b-a).toVector3D(), (p1-a).toVector3D());
	QVector3D cp2 = QVector3D::crossProduct((b-a).toVector3D(), (p2-a).toVector3D());
	if(QVector3D::dotProduct(cp1, cp2) >= 0) return false;
	return true;
}

QVector2D RotateCW90(QVector2D p){
	return QVector2D(p.y(), -p.x());
}

QVector2D RotateCCW90(QVector2D p){
	return QVector2D(-p.y(), p.x());
}


