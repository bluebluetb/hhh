#ifndef DYNAMICMESH_IMPL_H
#define DYNAMICMESH_IMPL_H

#include "dynamicmesh.h"
#include <QPair>
#include <QMap>
#include <QDebug>
#include "geom.h"

namespace DMesh{
class Vertex: public MeshVertex{
	public:
		Vertex():MeshVertex(){real=false;id=nextId;nextId++;}
		Vertex(double x, double y, bool _real=false) : MeshVertex(x,y){Backup();real=_real;id=nextId;nextId++;}
		Vertex(QVector2D v, bool _real=false): MeshVertex(v){real=_real;id=nextId;nextId++;}
		~Vertex(){
			Q_ASSERT(true);
		}
		void RenewID(){id=nextId;nextId++;}
		QVector2D TextureCoords() const{return QVector2D(tx,ty);}
		QVector2D GetOrigCoords() const{return (hasOrigCoords)?origCoords:QVector2D(0,0);}
		bool HasOrigCoords() const{return hasOrigCoords;}
		void SetOrigCoords() {origCoords=toQVector2D();hasOrigCoords=true;}

		QVector<int> triangles;
		double value=0;//real data may have some value associated with it
		float tx,ty;
		QVector2D toQVector2D()const {return QVector2D(x(),y());}
		bool IsReal()const{return real;}
		void SetReal(bool val){real=val;}
		int ID() const{return id;}
		bool MoveTo(double x, double y, double hardMinHeightFactor, bool checkonly=false){return MoveTo(QVector2D(x,y),hardMinHeightFactor, checkonly);}
		bool MoveTo(QVector2D newPos, double hardMinHeightFactor, bool checkOnly=false);
		void SetName(QString _name){
			name=_name.toStdString();
		}
		void SetName(const char* _name){name=_name;}
		const unsigned char* GetName(){return (const unsigned char*)name.c_str();}
		QVector2D desiredPos=QVector2D(x(),y());
		void Backup(){backup.setX(x());backup.setY(y());}
		void Revert(){setX(backup.x());setY(backup.y());}
	private:
		std::string name="";
		bool real; //represents a non-grid point
		bool hasOrigCoords=false;
		QVector2D origCoords;
		int id;
		static int nextId;
		inline bool SegmentLineIntersect(QVector2D p1, QVector2D p2, QVector2D a, QVector2D b){
			QVector3D cp1 = QVector3D::crossProduct((b-a).toVector3D(), (p1-a).toVector3D());
			QVector3D cp2 = QVector3D::crossProduct((b-a).toVector3D(), (p2-a).toVector3D());
			if(QVector3D::dotProduct(cp1, cp2) >= 0) return false;
			return true;
		}
		inline QVector2D ProjectPointLine ( QVector2D a , QVector2D b , QVector2D p ) {
			return a + (b - a ) *QVector2D::dotProduct(( p - a ), (b - a )) /QVector2D::dotProduct(( b - a ),(b - a ));
		}
		QVector2D backup;
};
}

class ArrayMapper{
	public:
		ArrayMapper(): which(0){}
		ArrayMapper(QVector<QVector2D> *p, unsigned int i): p(p),which(i){}
		ArrayMapper(QVector<Vertex*> *pv, unsigned int i): vertex(true), pv(pv),which(i){}

		void Set(const unsigned i, double val){
			if(vertex){
				(which==0)?(*pv)[i]->setX(val):(*pv)[i]->setY(val);
			}
			else{
				(which==0)?(*p)[i].setX(val):(*p)[i].setY(val);
			}

		}
		double operator()(const unsigned i) const{
			if(vertex){
				return (which==0)?(*pv)[i]->x():(*pv)[i]->y();
			}
			else{
				return (which==0)?(*p)[i].x():(*p)[i].y();
			}
		}

		int size() const{return (vertex)?pv->size():p->size();}
	private:
		bool vertex=false;
		QVector<QVector2D> *p=nullptr;
		QVector<Vertex*> *pv=nullptr;
		unsigned int which;
};

class Vector2DContainer{
	public:
		virtual ArrayMapper & operator[](const unsigned i)=0;
		virtual double operator()(int xy, int i) const=0;
		virtual void resize(int n)=0;
		virtual int size()const=0;
		virtual void Set(unsigned i, double a, double b){x.Set(i,a);y.Set(i,b);}
		virtual bool MoveTo(unsigned i, double a, double b)=0;
	protected:
		ArrayMapper x,y;
};

class QVector2DContainer : public Vector2DContainer{
	public:
		QVector2DContainer(){vec=new QVector<QVector2D>();x=ArrayMapper(vec,0);y=ArrayMapper(vec,1);}
		virtual ~QVector2DContainer(){delete vec;}
		ArrayMapper & operator[](const unsigned i){return (i==0)?x:y;}
		double operator()(int xy, int i) const{return (xy==0)?(*vec)[i].x():(*vec)[i].y();}
		void resize(int n){vec->resize(n);}
		int size()const {return vec->size();}
		void append(QVector2D v){vec->append(v);}
		bool MoveTo(unsigned i, double a, double b){Set(i,a,b);return true;}
	private:
		QVector<QVector2D> * vec;
};

class VertexContainer: public Vector2DContainer {
	public:
		VertexContainer(double hardMinHeightFactor):hardMinHeightFactor(hardMinHeightFactor){vec=new QVector<Vertex*>();x=ArrayMapper(vec,0);y=ArrayMapper(vec,1);}
		~VertexContainer(){delete vec;}
		ArrayMapper & operator[](const unsigned i){return (i==0)?x:y;}
		double operator()(int xy, int i) const{return (xy==0)?(*vec)[i]->x():(*vec)[i]->y();}
		void resize(int n){vec->resize(n);}
		int size()const {return vec->size();}
		void append(Vertex* v){vec->append(v);}
		bool IsReal(unsigned i)const {return (*vec)[i]->IsReal();}
		int ID(unsigned i)const {return (*vec)[i]->ID();}
		bool MoveTo(unsigned i, double a, double b){
			Vertex* v=(*vec)[i];
			return v->MoveTo(a,b,hardMinHeightFactor);
		}
	private:
		double hardMinHeightFactor=0;
		QVector<Vertex*> * vec;
};

namespace DMesh{
#define TNOCOLOR 0
#define TGREEN 1
#define TRED 2
#define TBLUE 3


inline QVector2D ProjectPointLine ( QVector2D a , QVector2D b , QVector2D p ) {
	return a + (b - a ) *QVector2D::dotProduct(( p - a ), (b - a )) /QVector2D::dotProduct(( b - a ),(b - a ));
}


class Triangle: public MeshTriangle{
	public:
		Triangle():MeshTriangle(0,0,0){}
		Triangle(Vertex *a, Vertex *b, Vertex *c):MeshTriangle(a,b,c){UpdateMinHeight();}
		void AddToQVector(QVector<QVector2D> &r) const;
		TriangleP toTriangleP() const;
		int color=TNOCOLOR; //mostly used for debugging purposes
		void Set(unsigned i,MeshVertex* v){MeshTriangle::Set(i,v);UpdateMinHeight();}
		double OrigMinHeight(int v=-1){
			if(!minHeightSet){UpdateMinHeight();}
			if(v==-1){return minMinHeight;}
			return minHeight[v];
		}
		int stressEdge=-1;
		bool changed=true;
	private:
		double minHeight[3];
		bool minHeightSet=false;
		double minMinHeight=0;
		void UpdateMinHeight(){
			for(int i=0;i<3;i++){
				if(!((Vertex*)(*this)(i))->HasOrigCoords()) return;
			}
			TriangleP origTriangle(((Vertex*)a)->GetOrigCoords(), ((Vertex*)b)->GetOrigCoords(), ((Vertex*)c)->GetOrigCoords());
			for(int i=0;i<3;i++){
				minHeight[i]=(origTriangle(i)-ProjectPointLine(origTriangle((i+1)%3),origTriangle((i+2)%3),origTriangle(i))).length();
				if(minHeight[i]<minMinHeight || i==0){minMinHeight=minHeight[i];}
			}
			minHeightSet=true;
		}
};

class DynamicMesh: public _DynamicMesh<Vertex,Triangle>{
	public:
		DynamicMesh():_DynamicMesh(){D=QMap<QPair<int,int>, QPair<double,double> >();}
		DynamicMesh(QVector<QVector2D> p, QMap<QPair<int,int>, QPair<double,double> > &_D):_DynamicMesh(p){D=_D;}
		Vertex * GetVertexPointer(){return V;}
		void EdgeAdded(Vertex *x, Vertex *y){
			double d;
			if(x->IsReal() && y->IsReal()){return;}
			if(x->HasOrigCoords() && y->HasOrigCoords()){
				d=(x->GetOrigCoords()-y->GetOrigCoords()).length();
			}else{
				d=((*x)-(*y)).length();
			}
			if(x->ID()<y->ID()){
				D[qMakePair(x->ID(),y->ID())].first=d;
			}
			else{
				D[qMakePair(y->ID(),x->ID())].first=d;
			}
		}
		void EdgeRemoved(Vertex *x, Vertex *y){ //not always called by base class atm!
			bool stillUsed=false;
			for(int i=0;i<x->Data().size();i++){
				if(x->Data()[i]->Contains(y)){stillUsed=true;}
			}
			if(!stillUsed){
				D.remove(qMakePair(x->ID(),y->ID()));
			}
		}

		void TriangleAdded(Triangle *){

		}

		double Dis(int a, int b){
			QPair<int,int> pair;
			if(a<b){pair=qMakePair(a,b);}
			else{pair=qMakePair(b,a);}
			auto it=D.find(pair);
			auto max=[](double a,double b){return (a>b)?a:b;};
			if(it==D.end()){return 0;}
			else{return max(0,it.value().first+it.value().second);}
		}
		void SetDis(int a, int b, double val){
			if(a<b){D[qMakePair(a,b)].first=val;}
			else{D[qMakePair(b,a)].first=val;}
		}
		void SetExtraDis(int a, int b, double val){
			if(a<b){D[qMakePair(a,b)].second=val;}
			else{D[qMakePair(b,a)].second=val;}
		}
		double ExtraW(int a, int b){
			QPair<int,int> pair;
			if(a<b){pair=qMakePair(a,b);}
			else{pair=qMakePair(b,a);}
			auto it=extraW.find(pair);
			if(it==extraW.end()){return 0;}
			else{return it.value();}
		}

		void SetExtraW(int a, int b, double val){
			if(a<b){extraW[qMakePair(a,b)]=val;}
			else{extraW[qMakePair(b,a)]=val;}
		}

		void AddVertex(Vertex *v){
			_DynamicMesh<Vertex,Triangle>::AddVertex(v);
		}
		virtual void Triangulate(){
			_DynamicMesh<Vertex,Triangle>::Triangulate();
			Vertex * itVert=V;
			while(itVert!=nullptr){
				itVert->SetOrigCoords();
				itVert=(Vertex*)itVert->next;
			}
		}
		QMap<QPair<int,int>,QPair<double,double> > D;
		QMap<QPair<int,int>,double> extraW;
};

typedef _TriangleIterator<Vertex,Triangle> TriangleIterator;
typedef _VertexIterator<Vertex,Triangle> VertexIterator;
}
#endif // DYNAMICMESH_IMPL_H
