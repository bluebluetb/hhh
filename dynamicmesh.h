/* DynamicMesh class (Work in progress)
 * Provides traversal, edge flipping, triangle subdivision and triangle merging
 * Check bottom of file for example implementation file with correct inheritance
 *
 * The vertices in a triangle are stored in clockwise order.
 *
 * Using a c++ interface to the Triangle package for Delaunay triangulations.
 * The Triangle package was written by Jonathan Richard Shewchuk. (Details: triangle.h)
 * The C++ interface by Piyush Kumar and slightly extended by me. (Details: del_interface.h)
 *
 * Author: Quirijn Bouts - www.qbouts.com
 */

#ifndef DYNAMICMESH_H
#define DYNAMICMESH_H

#include <QVector>
#include <QVector2D>
#include <QVector3D>
#include <QPair>
#include <QMap>
#include <QPolygonF>

#include <limits>

#include "del_interface.h"

namespace DMesh{

template<class T>
class UnorderedTriple{
	public:
		//constructors
		UnorderedTriple(): a(0),b(0),c(0){}
		UnorderedTriple(T a, T b, T c): a(a),b(b),c(c){}
		//accessors
		T operator()(const unsigned i) const{return (i==0)?a:(i==1)?b:c;}
		//functions
		bool HasDuplicates()const {return !(a!=b && a!=c && c!=b);}
		bool Contains(T val) const{for(int i=0;i<3;i++){if((*this)(i)==val){return true;}}return false;}
		bool Equals(const UnorderedTriple &t2) const{return (t2.Contains((*this)(0)) && t2.Contains((*this)(1)) && t2.Contains((*this)(2)));}
		QVector<T> ToQVector(){QVector<T> v; AddToVector(v); return v;}
		void AddToVector(QVector<T> &vec){vec.append(a); vec.append(b); vec.append(c); return;}
	protected:
		virtual T &operator[](const unsigned i) {return (i==0)?a:(i==1)?b:c;}
		T a,b,c;
};

template<class DATA>
class DataVertex : public QVector2D{
	public:
		//constructors
		DataVertex():QVector2D(0,0){}
		DataVertex(double x, double y): QVector2D(x,y){}
		DataVertex(double x, double y, DATA &data):QVector2D(x,y),data(data){}
		DATA &Data(){return data;}
	private:
		DATA data;
};

class MeshTriangle;

class MeshVertex : public DataVertex<QVector<MeshTriangle*> >{
	public:
		MeshVertex():DataVertex<QVector<MeshTriangle*> >(){}
		MeshVertex(QVector2D v): DataVertex<QVector<MeshTriangle*> >(v.x(),v.y()){}
		MeshVertex(double x, double y): DataVertex<QVector<MeshTriangle*> >(x,y){}
		MeshVertex(QVector2D v, QVector<MeshTriangle*> &data):DataVertex<QVector<MeshTriangle*> >(v.x(),v.y(),data){}
		MeshVertex(double x, double y, QVector<MeshTriangle*> &data):DataVertex<QVector<MeshTriangle*> >(x,y,data){}
		double crossProduct(const QVector2D &other){return (x()*other.y()-y()*other.x());} //lower camelcase to match QVector::dotProduct
		MeshVertex* next=nullptr;
		MeshVertex* prev=nullptr;
		bool RemoveFromTriangle(MeshTriangle* T);
		bool IsInTriangle(MeshTriangle* T){return Data().contains(T);}
};

class Edge: public QPair<QVector2D,QVector2D>{
	public:
		Edge():QPair(QVector2D(0,0),QVector2D(0,0)){}
		Edge(QVector2D a, QVector2D b){first=a; second=b;}
		double Length(){return (first-second).length();}
		double LengthSquared(){return (first-second).lengthSquared();}
};

class TriangleP: public UnorderedTriple<QVector2D>{
	public:
		TriangleP():UnorderedTriple(QVector2D(0,0),QVector2D(0,0),QVector2D(0,0)){}
		TriangleP(QVector2D a, QVector2D b, QVector2D c):UnorderedTriple(a,b,c){}
		Edge GetEdge(const unsigned i) const{return Edge((*this)(i),(*this)((i+1)%3));}
		bool Encloses(const QVector2D &p) const;
		double MinHeight() const;
		double Area() const;
		bool IsClockwise() const;
		QVector2D Centroid() const{return ((*this)(0)+(*this)(1)+(*this)(2))/3;}
		void Scale(double factor);
		QVector2D &operator[](const unsigned i) {return (i==0)?a:(i==1)?b:c;}
};

class MeshTriangle: public UnorderedTriple<MeshVertex*>{
	public:
		MeshTriangle(MeshVertex* a, MeshVertex* b, MeshVertex* c):UnorderedTriple<MeshVertex*>(a,b,c){Q_ASSERT(a!=nullptr && b!=nullptr && c!=nullptr);}
		MeshTriangle* next=nullptr;
		MeshTriangle* prev=nullptr;
		virtual void Set(unsigned int i, MeshVertex* v);
		virtual bool Swap(MeshVertex * oldV, MeshVertex * newV);
		TriangleP SwapResult(MeshVertex * oldV, MeshVertex * newV) const; //shows the result of a swap but does not actually modify the meshtriangle
		int GetPos(MeshVertex* v)const {int pos=-1;for(int i=0;i<3;i++){if((*this)(i)==v) pos=i;}return pos;}
};

template<class VERTEX,class TRIANGLE>
class _DynamicMesh; //forward declaration

template<class VERTEX,class TRIANGLE>
class _VertexIterator{
	public:
		//constructors
		_VertexIterator(): mesh(nullptr),V(nullptr){}
		_VertexIterator(_DynamicMesh<VERTEX,TRIANGLE> * mesh, VERTEX * V): mesh(mesh),V(V){}
		//equality tests
		bool operator ==(const _VertexIterator &other)const {return (V==other.V);}
		bool operator !=(const _VertexIterator &other)const {return (V!=other.V);}
		//const accessor
		inline double x()const {return V->x();}
		inline double y()const {return V->y();}
		QVector2D GetV() const{return QVector2D(V->x(), V->y());}
		//accessors
		VERTEX &Vertex(){return *V;}
		VERTEX *VertexPointer(){return V;}
		//traversal
		void Next(){if(V==End()){return;}V=(VERTEX*)V->next;}
		void Prev(){if(V==End()){return;}V=(VERTEX*)V->prev;}
		//functions
		bool Last() const{return (V!=nullptr && V->next==End());}
		bool First() const{return (V!=nullptr && V->prev==End());}
		int nNeighbours() const{return V->Data().size();}
		void GetNeighbours(QVector<VERTEX*> &N) const;
		//mesh modifications
		void RemoveAndMerge(); //removes vertex and retriangulates the resulting gap, afterwards the iterator will point to the next vertex
	private:
		_DynamicMesh<VERTEX,TRIANGLE> * mesh;
		VERTEX * V;
		inline VERTEX * End() const {return nullptr;}
};


template<class VERTEX,class TRIANGLE>
class _TriangleIterator{
	public:
		//constructors
		_TriangleIterator();
		_TriangleIterator(_DynamicMesh<VERTEX,TRIANGLE> * mesh, TRIANGLE * T);
		//equality tests
		bool operator ==(const _TriangleIterator &other)const;
		bool operator !=(const _TriangleIterator &other)const;
		//const accessors
		TRIANGLE &GetTriangle()const;
		TriangleP GetTriangleP()const;
		VERTEX &GetVertexA()const;
		VERTEX &GetVertexB()const;
		VERTEX &GetVertexC()const;
		VERTEX &operator()(const unsigned i) const;
		Edge GetEdge(const unsigned i) const;
		//accessors
		TRIANGLE &Triangle();
		VERTEX &operator[](const unsigned i);
		_TriangleIterator<VERTEX,TRIANGLE> GetNeighbourIt(const unsigned a, const unsigned b);
		_TriangleIterator<VERTEX,TRIANGLE> GetNeighbourIt(const unsigned e);
		//traversal
		void Next();
		void Prev();
		void Neighbour(const unsigned a, const unsigned b);
		void Neighbour(const unsigned e);
		//some convenience functions
		void NeighbourAB();
		void NeighbourBA();
		void NeighbourBC();
		void NeighbourCB();
		void NeighbourAC();
		void NeighbourCA();

		//functions
		bool Last() const{return (T!=nullptr && T->next==End());}
		bool First() const{return (T!=nullptr && T->prev==End());}
		bool HasNeighbour(const unsigned a, const unsigned b);
		bool HasNeighbour(const unsigned e);
		//some convenience functions
		bool HasNeighbourAB();
		bool HasNeighbourBA();
		bool HasNeighbourBC();
		bool HasNeighbourCB();
		bool HasNeighbourAC();
		bool HasNeighbourCA();

		//mesh modifications
		bool FlipEdge(const unsigned i); //return true unless no triangle neighbours edge i
		template<typename SUBDIVFUNCTION>
		int SubDivide(const unsigned e, SUBDIVFUNCTION F); //returns number of subdivisions (1 or 2 if edde is shared by two triangles)
		int SubDivide(const unsigned e){return SubDivide(e,[](VERTEX * a, VERTEX * b, TRIANGLE*){return new VERTEX(((*a)+(*b))/2);});}
		void Merge(const unsigned v); //removes the vertex and retriangulates the resulting gap, afterwards the iterator will point to the next triangle which was also in the original mesh
	private:
		_DynamicMesh<VERTEX,TRIANGLE> * mesh;
		TRIANGLE * T;
		inline TRIANGLE * End() const;
		//functions
		TRIANGLE * NeighbourPointer(const unsigned e);
		TRIANGLE * NeighbourPointer(const unsigned a, const unsigned b);
};

template<class VERTEX,class TRIANGLE>
class _DynamicMesh
{
	public:
		//constructors
		_DynamicMesh();
		_DynamicMesh(const QVector<QVector2D> &vertices);
		//copy constructor
		_DynamicMesh(const _DynamicMesh<VERTEX,TRIANGLE> & other);
		//assignment
		void operator=(const _DynamicMesh<VERTEX,TRIANGLE> & other);
		//descructor
		virtual ~_DynamicMesh();
		//accessors
		_TriangleIterator<VERTEX,TRIANGLE> Begin(){return _TriangleIterator<VERTEX,TRIANGLE>(this,T);}
		_TriangleIterator<VERTEX,TRIANGLE> End(){return _TriangleIterator<VERTEX,TRIANGLE>(this,nullptr);}
		_TriangleIterator<VERTEX,TRIANGLE> Iterator(TRIANGLE * t) {return _TriangleIterator<VERTEX,TRIANGLE>(this,t);}
		_VertexIterator<VERTEX,TRIANGLE> VBegin(){return _VertexIterator<VERTEX,TRIANGLE>(this,V);}
		_VertexIterator<VERTEX,TRIANGLE> VEnd(){return _VertexIterator<VERTEX,TRIANGLE>(this,nullptr);}
		_VertexIterator<VERTEX,TRIANGLE> VIterator(VERTEX * v) {return _VertexIterator<VERTEX,TRIANGLE>(this,v);}
		VERTEX GetVertex(const unsigned i)const {return V[i];}
		TRIANGLE GetTriangle(const unsigned i)const {return T[i];}
		//functions
		int nPoints()const{return vertexCount;}
		int nTriangles()const{return triangleCount;}
		virtual void SetVertices(const QVector<QVector2D> &vertices);//also triggers a Triangulate() call
		virtual void SetVertices(const QVector<VERTEX> &vertices);//also triggers a Triangulate() call
		virtual void SetVertices(VERTEX * vertices);//also triggers a Triangulate() call

		//virtual functions
		virtual void EdgeAdded(VERTEX*,VERTEX*){} //called when an edge is added
		virtual void EdgeRemoved(VERTEX*,VERTEX*){} //called when an edge is removed through a flip only!!
		virtual void TriangleAdded(TRIANGLE *){}
		virtual void VertexAdded(VERTEX *){}

	protected:
		friend class _TriangleIterator<VERTEX,TRIANGLE>;
		friend class _VertexIterator<VERTEX,TRIANGLE>;
		VERTEX * V=nullptr;
		VERTEX * lastV=nullptr;
		TRIANGLE * T=nullptr;
		TRIANGLE * lastT=nullptr;
		int triangleCount=0,vertexCount=0;
		virtual void AddVertex(VERTEX *v);
		virtual void DeleteVertex(VERTEX *v);
		virtual TRIANGLE* AddTriangle(VERTEX * a, VERTEX * b, VERTEX * c); //adds triangle at also adds its index to the vertices involved
		virtual void DeleteTriangle(TRIANGLE* t, bool deleteDisconnectedVertices=true); //removes triangle and removes its index from vertices involved
		virtual void Triangulate(); //performs a delaunay triangulation on V
		void CopyValues(const _DynamicMesh<VERTEX,TRIANGLE> & other);
};
}

/*************************************************/
/**************** IMPLEMENTATION *****************/
/*************************************************/

using namespace DMesh;

namespace DMesh{
inline bool SameSide(const QVector2D &p1, const QVector2D &p2,const QVector2D &a,const QVector2D &b){
	QVector3D cp1 = QVector3D::crossProduct((b-a).toVector3D(), (p1-a).toVector3D());
	QVector3D cp2 = QVector3D::crossProduct((b-a).toVector3D(), (p2-a).toVector3D());
	if(QVector3D::dotProduct(cp1, cp2) >= 0) return true;
	return false;
}
}

/*****   Vertex Iterator    *****/
template<class VERTEX,class TRIANGLE>
void _VertexIterator<VERTEX, TRIANGLE>::GetNeighbours(QVector<VERTEX *> &N) const{
	for(int i=0;i<V->Data().size();i++){
		int j=V->Data()[i]->GetPos(V);
		Q_ASSERT(j!=-1);
		VERTEX* v=(VERTEX*)(*(V->Data()[i]))((j+1)%3);
		if(!N.contains(v)){N.append(v);}
		v=(VERTEX*)(*(V->Data()[i]))((j+2)%3);
		if(!N.contains(v)){N.append(v);}
	}
	return;
}

/*****   TriangleIterator   *****/
//constructors
template<class VERTEX,class TRIANGLE>
_TriangleIterator<VERTEX,TRIANGLE>::_TriangleIterator():mesh(nullptr), T(nullptr){}
template<class VERTEX,class TRIANGLE>
_TriangleIterator<VERTEX,TRIANGLE>::_TriangleIterator(_DynamicMesh<VERTEX,TRIANGLE> *mesh, TRIANGLE *T): mesh(mesh), T(T){}

//equality tests
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::operator ==(const _TriangleIterator &other)const{
	return (T==other.T);
}

template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::operator !=(const _TriangleIterator &other)const{
	return (T!=other.T);
}

//const accessors
template<class VERTEX,class TRIANGLE>
TRIANGLE &_TriangleIterator<VERTEX,TRIANGLE>::GetTriangle()const {return *(TRIANGLE*)(T);}
template<class VERTEX,class TRIANGLE>
TriangleP _TriangleIterator<VERTEX,TRIANGLE>::GetTriangleP()const {return TriangleP(*(*T)(0),*(*T)(1),*(*T)(2));}
template<class VERTEX,class TRIANGLE>
VERTEX &_TriangleIterator<VERTEX, TRIANGLE>::GetVertexA()const {return (*(*T)[0]);}
template<class VERTEX,class TRIANGLE>
VERTEX &_TriangleIterator<VERTEX,TRIANGLE>::GetVertexB()const {return (*(*T)[1]);}
template<class VERTEX,class TRIANGLE>
VERTEX &_TriangleIterator<VERTEX,TRIANGLE>::GetVertexC()const {return (*(*T)[2]);}
template<class VERTEX,class TRIANGLE>
VERTEX &_TriangleIterator<VERTEX,TRIANGLE>::operator()(const unsigned i) const{return (*(VERTEX*)((*T)(i)));}
template<class VERTEX,class TRIANGLE>
Edge _TriangleIterator<VERTEX,TRIANGLE>::GetEdge(const unsigned i) const{return GetTriangleP().GetEdge(i);}

//accessors
template<class VERTEX,class TRIANGLE>
TRIANGLE & _TriangleIterator<VERTEX,TRIANGLE>::Triangle(){return *T;}
template<class VERTEX,class TRIANGLE>
VERTEX &_TriangleIterator<VERTEX,TRIANGLE>::operator[](const unsigned i){return (*((VERTEX*)(*T)(i)));}
template<class VERTEX,class TRIANGLE>
_TriangleIterator<VERTEX,TRIANGLE> _TriangleIterator<VERTEX,TRIANGLE>::GetNeighbourIt(const unsigned a, const unsigned b){
	return _TriangleIterator<VERTEX,TRIANGLE>(mesh,NeighbourPointer(a,b));
}
template<class VERTEX,class TRIANGLE>
_TriangleIterator<VERTEX,TRIANGLE> _TriangleIterator<VERTEX,TRIANGLE>::GetNeighbourIt(const unsigned e){
	return _TriangleIterator<VERTEX,TRIANGLE>(mesh,NeighbourPointer(e));
}
//traversal
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::Next(){
	if(T==End()){return;}
	T=(TRIANGLE*)(T->next);
}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::Prev(){
	if(T==End()){return;}
	T=(TRIANGLE*)(T->prev);
}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::Neighbour(const unsigned a, const unsigned b){
	T=NeighbourPointer(a,b);
}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::Neighbour(const unsigned e){
	T=NeighbourPointer(e);
}
//some convenience functions
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::NeighbourAB(){Neighbour(0,1);}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::NeighbourBA(){Neighbour(0,1);}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::NeighbourBC(){Neighbour(1,2);}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::NeighbourCB(){Neighbour(1,2);}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::NeighbourAC(){Neighbour(0,2);}
template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX,TRIANGLE>::NeighbourCA(){Neighbour(0,2);}

//functions
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbour(const unsigned a, const unsigned b){
	return (NeighbourPointer(a,b)!=End());
}
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbour(const unsigned e){
	return (NeighbourPointer(e)!=End());
}
//some convenience functions
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbourAB(){return HasNeighbour(0,1);}
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbourBA(){return HasNeighbour(0,1);}
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbourBC(){return HasNeighbour(1,2);}
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbourCB(){return HasNeighbour(1,2);}
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbourAC(){return HasNeighbour(0,2);}
template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::HasNeighbourCA(){return HasNeighbour(0,2);}


//private functions
template<class VERTEX,class TRIANGLE>
TRIANGLE * _TriangleIterator<VERTEX,TRIANGLE>::End() const{return nullptr;}

template<class VERTEX,class TRIANGLE>
TRIANGLE *_TriangleIterator<VERTEX, TRIANGLE>::NeighbourPointer(const unsigned e){
	return NeighbourPointer(e, (e+1)%3);
}

template<class VERTEX,class TRIANGLE>
TRIANGLE *_TriangleIterator<VERTEX, TRIANGLE>::NeighbourPointer(const unsigned a, const unsigned b){
	Q_ASSERT(a<=3 && b<=3);
	//for all triangles sharing vertex a, only the correct neighbour also contains b
	for(int i=0;i<(*T)(a)->Data().size();i++){
		MeshTriangle * t2=(*T)(a)->Data()[i];
		if(t2==T){continue;}//don't compare with yourself
		if(t2->Contains((*T)(b))){return (TRIANGLE*)t2;break;}
	}
	//no neighbour sharing edge
	return End();
}

template<class T>
int RemoveValue(QVector<T> &vec, const T &value){ //removes first occurence of value from vec, returns index, or -1 if no occurence
	for(int i=0;i<vec.size();i++){
		if(vec[i]==value){vec.remove(i);return i;}
	}
	return -1;
}

template<class VERTEX,class TRIANGLE>
bool _TriangleIterator<VERTEX,TRIANGLE>::FlipEdge(const unsigned e){ //return true unless no triangle neighbours edge i / resulting quad is concave
	if(!HasNeighbour(e)){return false;}
	TRIANGLE &one = Triangle();
	_TriangleIterator<VERTEX,TRIANGLE> it2=GetNeighbourIt(e);
	TRIANGLE &two = it2.Triangle();
	QPair<int,int> e2; //indices of the flip edge in triangle two, first corresponds to e, second to (e+1)%3
	int otherp=-1; //the other point of triangle two
	for(int i=0;i<3;i++){
		if(two(i)==one(e) && two((i+1)%3)==one((e+1)%3)){
			e2=qMakePair(i,(i+1)%3);otherp=(i+2)%3;break;
		}
		if(two(i)==one((e+1)%3) && two((i+1)%3)==one(e)){
			e2=qMakePair((i+1)%3,i);otherp=(i+2)%3;break;
		}
	}
	Q_ASSERT(otherp != -1);
	//check if quad is convex
	QVector2D p[4];
	p[0]=(*one((e+1)%3));
	p[1]=(*one((e+2)%3));
	p[2]=(*one((e)%3));
	p[3]=(*two(otherp));
	bool det=(VERTEX(p[3]-p[0]).crossProduct(p[1]-p[0])>0);
	for(int i=1;i<4;i++){
		if((VERTEX(p[(i+3)%4]-p[i]).crossProduct(p[(i+1)%4]-p[i])>0)!=det){return false;}
	}
	//its safe to flip!
	Q_ASSERT(one(e)==two(e2.first) && one((e+1)%3)==two(e2.second));
	//fix vertex data
	/*if(RemoveValue(one[e]->Data(),(MeshTriangle*)T)==-1){Q_ASSERT(false);}
	if(RemoveValue(one[(e+1)%3]->Data(),(MeshTriangle*)it2.T)==-1){Q_ASSERT(false);}
	one[(e+2)%3]->Data().append((MeshTriangle*)it2.T);
	two[otherp]->Data().append((MeshTriangle*)T);*/
	//do the actual flip
	one.Set(e,two(otherp));
	two.Set(e2.second,one((e+2)%3));
	//call edge modification functions
	mesh->EdgeRemoved((VERTEX*)two(e2.first),(VERTEX*)one((e+1)%3));
	mesh->EdgeAdded((VERTEX*)one((e+2)%3),(VERTEX*)two(otherp));//*/
	return true;
}
/*   /|\              / \
 *	/4|2\            / 2 \
 * e--v--e+1  <---  e-----e+1
 *  \3|1/            \ 1 /
 *   \|/              \ /
 *
 */
template<class VERTEX,class TRIANGLE> template<typename SUBDIVFUNCTION>
int _TriangleIterator<VERTEX,TRIANGLE>::SubDivide(const unsigned e, SUBDIVFUNCTION F){
	TRIANGLE &one=Triangle();
	_TriangleIterator<VERTEX,TRIANGLE> it2=GetNeighbourIt(e);
	bool hasNeighbour=(it2!=mesh->End());
	VERTEX * v=F((VERTEX*)one(e),(VERTEX*)one((e+1)%3),T);
	mesh->AddVertex(v);
	//RemoveValue(one(e)->Data(),(MeshTriangle*)T);
	//v->Data().append(T);
	mesh->AddTriangle((VERTEX*)one((e+2)%3),(VERTEX*)one(e),v); //also adds it to the vertex data

	VERTEX * oldOne_e=(VERTEX*)one(e);
	one.Set(e,v);
	//finished subdividing this side, check if there is a neighbour to fix
	if(hasNeighbour){
		VERTEX * otherp=nullptr;
		//v->Data().append(it2.T);
		for(int i=0;i<3;i++){
			if(it2.Triangle()(i)!=oldOne_e && it2.Triangle()(i)!=one((e+1)%3)){
				otherp=(VERTEX*)it2.Triangle()(i);
			}
		}
		for(int i=0;i<3;i++){
			if(it2.Triangle()(i)==oldOne_e){
				//RemoveValue(oldOne_e->Data(),(MeshTriangle*)it2.T);
				it2.Triangle().Set(i,v);
				break;
			}
		}
		Q_ASSERT(otherp!=nullptr);
		//and add the new triangle
		mesh->AddTriangle(otherp,v,oldOne_e);
		return 2;
	}
	return 1;
}

/*     /|\             / \
 *    /3|4\           / 3 \
 * e+2--e---   --->  -------
 *    \1|2/           \ 1 /
 *     \|/             \ /
 *     e+1
 */

inline bool IsLeft(QPointF a, QPointF b, QPointF c){
	return ((b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x())) > 0;
}

template<class VERTEX, class TRIANGLE>
void _VertexIterator<VERTEX,TRIANGLE>::RemoveAndMerge(){
	QVector<VERTEX*> verts;
	GetNeighbours(verts);
	VERTEX * rm=this->V;
	qSort(verts.begin(),verts.end(),[rm](VERTEX* a, VERTEX* b){
		QVector2D aa=(*a)-(*rm);
		QVector2D bb=(*b)-(*rm);
		return atan2(aa.y(),aa.x())<atan2(bb.y(),bb.x());
	});
	//by sorting on angle around rm, we know that there is an edge between any verts[i],verts[i+1] pair

	for(int i=V->Data().size()-1;V->Data().size()>0;i=V->Data().size()-1){
		mesh->DeleteTriangle((TRIANGLE*)(V->Data()[i]),false);
	}
	mesh->DeleteVertex(rm);
	if(verts.size()==2){return;}
	std::vector<Delaunay::Point> points;
	QVector<VERTEX*> vPointers;
	std::vector<Delaunay::Segment> segments;
	for(int i=0;i<verts.size();i++){
		points.push_back(Delaunay::Point(verts[i]->x(),verts[i]->y()));
		vPointers.push_back(verts[i]);
		segments.push_back(Delaunay::Segment(i,(i+1)%verts.size(),2));
	}
	Delaunay delobject(points,segments);
	delobject.Triangulate();

	if(delobject.nvertices()==points.size()){
		for (Delaunay::fIterator fit = delobject.fbegin(); fit != delobject.fend(); ++fit) {
			mesh->AddTriangle(vPointers[delobject.Dest(fit)],vPointers[delobject.Org(fit)], vPointers[delobject.Apex(fit)]);
		}
	}
	else{ //in some rare cases Triangle adds an extra point, since this causes problems for us, we use the ear clipping method te get a different triangulation
		//uses earClipping method (not optimal, but it does the job)
		QVector<unsigned int> pi;//stores indiced of our clipped polygon
		QVector<QPointF> pts;
		for(unsigned int i=0;i<points.size();i++){pts.append(QPointF(points[i][0],points[i][1]));}
		for(int i=pts.size()-1;i>=0;i--){pi.push_back(i);}//we start with the complete polygon
		while(pi.size()>3){//while we have more than a triangle left
			//find pair of segments that has included angle < 180
			for(int i=0;i<pi.size()-2;i++){
				//check if point i+2 is to the "right" of line through i, i+1 (because of CW order this implies included angle < 180)
				if(!IsLeft(pts[pi[i]],pts[pi[i+1]],pts[pi[i+2]])){
					//now check if segment i, i+2 lies completely in the polygon (or if no other polygon points lie in triangle i,i+1,i+2)
					bool ear=true;
					QPolygonF tmp;
					tmp.append(pts[pi[i]]);tmp.append(pts[pi[i+1]]);tmp.append(pts[pi[i+2]]);
					for(int j=0;j<pts.size();j++){//j<Size
						if(j!=pi[i] && j!=pi[i+1] && j!=pi[i+2]){
							if(tmp.containsPoint(pts[j],Qt::OddEvenFill)){
								ear=false; //not a valid ear
								break;
							}
						}
					}
					if(ear){
						//i, i+1, i+2  is a valid ear and can be cut of
						mesh->AddTriangle(vPointers[pi[i]],vPointers[pi[i+1]],vPointers[pi[i+2]]); // add to triangulation
						pi.remove(i+1); //update polygon
						break; //time to start looking for the next ear
					}
				}
			}
		}
		//add last triangle
		mesh->AddTriangle(vPointers[pi[0]],vPointers[pi[1]],vPointers[pi[2]]);
	}
}

template<class VERTEX,class TRIANGLE>
void _TriangleIterator<VERTEX, TRIANGLE>::Merge(const unsigned v){
	auto it= _VertexIterator<VERTEX,TRIANGLE>(this->mesh,&(*this)[v]);
	it.RemoveAndMerge();

	/*Q_ASSERT(verts.size()>=2);
	if(verts.size()==2){
		//only delete triangle
		mesh->DeleteTriangle((TRIANGLE*)vertex->Data()[0]);//will also cleanup disconnected vertices
		return 0;
	}
	if(verts.size()==3){
		Q_ASSERT(vertex->Data().size()==2 || vertex->Data().size()==3);
		//easy case, only one triangle
		//delete other triangle(s) and expand the this one
		while(vertex->Data().size()>1){
			for(int i=0;i<vertex->Data().size();i++){
				TRIANGLE * t=(TRIANGLE*)vertex->Data()[i];
				if(t!=T){
					mesh->DeleteTriangle(t,false);
				}
			}
		}
		//T contains atleast 2 of the 3 vertices of the expanded triangle, find the other one
		VERTEX* newVertex=nullptr;
		for(int i=0;i<verts.size();i++){
			if(verts[i]!=(*T)((e+1)%3) && verts[i]!=(*T)((e+2)%3)){
				newVertex=verts[i];
			}
		}
		Q_ASSERT(newVertex!=nullptr);
		T->Swap((*T)(e),newVertex);
		//cleanup now disconnected vertex
		mesh->DeleteVertex(vertex);
		return 1;
	}
	if(verts.size()==4){
		//return -1;
		//find correct edge to insert
		for(int i=0;i<verts.size();i++){
			for(int j=i+1;j<verts.size();j++){
				int other1,other2;
				for(int k=0;k<verts.size();k++){
					if(k!=i && k!=j){other1=k;}
				}
				for(int k=0;k<verts.size();k++){
					if(k!=i && k!=j && k!=other1){other2=k;}
				}
				if(SameSide(*(verts[other1]),*(verts[other2]),*(verts[i]),*(verts[j]))){
					continue;
				}
				TriangleP one(*verts[i],*verts[j],*verts[other1]);
				TriangleP two(*verts[i],*verts[j],*verts[other2]);
				if(!one.Encloses(*verts[other2]) && !two.Encloses(*verts[other1])){
					//we found a good edge!
					while(vertex->Data().size()>0){
						mesh->DeleteTriangle((TRIANGLE*)vertex->Data()[0],false);
					}
					if(one.IsClockwise()){
						T=mesh->AddTriangle(verts[i],verts[j],verts[other1]);
					}else{
						T=mesh->AddTriangle(verts[j],verts[i],verts[other1]);
					}
					TRIANGLE *T2;
					if(two.IsClockwise()){
						T2=mesh->AddTriangle(verts[i],verts[j],verts[other2]);
					}else{
						T2=mesh->AddTriangle(verts[j],verts[i],verts[other2]);
					}
					mesh->DeleteVertex(vertex);
					return 2;
				}
			}
		}
	}
	if(verts.size()>4){


	}
	return -1;
	_TriangleIterator<VERTEX,TRIANGLE> it2=GetNeighbourIt(e);
	TRIANGLE &one=Triangle();
	Q_ASSERT(GetTriangleP().IsClockwise());
	Q_ASSERT(it2.GetTriangleP().IsClockwise());
	int eTriangles=one(e)->Data().size();
	if(it2==mesh->End()){
		return 0;
	}
	if(eTriangles!=2 && eTriangles!=4){
		return 0;
	}

	TRIANGLE &two=it2.Triangle();

	//info for first merge
	VERTEX* otherp=nullptr;
	for(int i=0;i<3;i++){
		if(two(i)!=one(e) && two(i)!=one((e+1)%3)){
			otherp=(VERTEX*)two(i);
		}
	}
	Q_ASSERT(otherp!=nullptr);

	_TriangleIterator<VERTEX,TRIANGLE> it3;
	_TriangleIterator<VERTEX,TRIANGLE> it4;
	if(eTriangles==4){
		for(int i=0;i<eTriangles;i++){
			TRIANGLE * t=(TRIANGLE*)one(e)->Data()[i];
			if(t!=T && t!=it2.T && !t->Contains(otherp)){
				it3=_TriangleIterator<VERTEX,TRIANGLE>(mesh,t);
			}
		}
		for(int i=0;i<eTriangles;i++){
			TRIANGLE * t=(TRIANGLE*)one(e)->Data()[i];
			if(t!=T && t!=it2.T && t!=it3.T){
				it4=_TriangleIterator<VERTEX,TRIANGLE>(mesh,t);
			}
		}
		Q_ASSERT(it3!=mesh->End());
		Q_ASSERT(it4!=mesh->End());
		Q_ASSERT(it2!=it4);
		Q_ASSERT(it3.GetTriangleP().IsClockwise());
		Q_ASSERT(it4.GetTriangleP().IsClockwise());

		//check if we can expand three without making it counterclockwise
		TriangleP newT3=it3.T->SwapResult(one(e),otherp);
		if(!newT3.IsClockwise()){
			//abort
			return 0;
			//this can happen if the four triangles form a concave quad
			//in such a case, merging along the wrong edge results in triangle flips
		}
		//expand three
		TriangleP t1=GetTriangleP();
		TriangleP t2=it2.GetTriangleP();
		TriangleP t3=it3.GetTriangleP();
		TriangleP t4=it4.GetTriangleP();
		bool ret=it3.T->Swap(one(e),otherp);
		Q_ASSERT(ret);
		TriangleP t32=it3.GetTriangleP();
		mesh->DeleteTriangle(it4.T);
		TriangleP t33=it3.GetTriangleP();
		Q_ASSERT(it3.GetTriangleP().IsClockwise());
	}
	//expand one
	T->Set(e,otherp);
	mesh->DeleteTriangle(it2.T); //should also disconnect the vertex and delete it
	Q_ASSERT(GetTriangleP().IsClockwise());
	if(eTriangles==4){return 2;}
	return 1;*/
}

/*****   DynamicMesh   *****/
template<class VERTEX,class TRIANGLE>
_DynamicMesh<VERTEX,TRIANGLE>::_DynamicMesh(){
}

template<class VERTEX,class TRIANGLE>
_DynamicMesh<VERTEX,TRIANGLE>::_DynamicMesh(const QVector<QVector2D> &vertices){
	SetVertices(vertices);
}

template<class VERTEX,class TRIANGLE>
_DynamicMesh<VERTEX,TRIANGLE>::_DynamicMesh(const _DynamicMesh<VERTEX,TRIANGLE> & other){
	CopyValues(other);
}

template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::operator=(const _DynamicMesh<VERTEX,TRIANGLE> & other){
	CopyValues(other);
}

template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::CopyValues(const _DynamicMesh<VERTEX,TRIANGLE> & other){
	//cleanup
	while(T!=nullptr){DeleteTriangle(T);}
	while(V!=nullptr){DeleteVertex(V);}
	//and assign
	QMap<VERTEX*,VERTEX*> oldNewMap; //map from original to copy
	MeshVertex * v=other.V;
	while(v!=nullptr){
		VERTEX* newV=new VERTEX(QVector2D(v->x(),v->y()));
		oldNewMap[(VERTEX*)v]=newV;
		AddVertex(newV);
		v=v->next;
	}
	MeshTriangle * t=other.T;
	typename QMap<VERTEX*,VERTEX*>::Iterator it;
	while(t!=nullptr){
		VERTEX *verts[3];
		for(int i=0;i<3;i++){
			it=oldNewMap.find((VERTEX*)(*t)(i));
			Q_ASSERT(it!=oldNewMap.end());
			verts[i]=it.value();
		}
		AddTriangle(verts[0],verts[1],verts[2]);
		t=t->next;
	}

}

template<class VERTEX,class TRIANGLE>
_DynamicMesh<VERTEX,TRIANGLE>::~_DynamicMesh(){
	while(T!=nullptr){DeleteTriangle(T);}
	while(V!=nullptr){DeleteVertex(V);}
}

template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::SetVertices(VERTEX * v){
	for(int i=0;i<vertexCount;i++){
		DeleteVertex(V);
	}
	Q_ASSERT(vertexCount==0);
	while(v!=nullptr){
		AddVertex(v);
		v=(VERTEX*)v->next;
	}
	Triangulate();
}

template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::SetVertices(const QVector<QVector2D> &vertices){
	for(int i=0;i<vertexCount;i++){
		DeleteVertex(V);
	}
	Q_ASSERT(vertexCount==0);

	for(int i=0;i<vertices.size();i++){
		VERTEX * v=new VERTEX(vertices[i]);
		AddVertex(v);
	}
	Triangulate();
}

template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::SetVertices(const QVector<VERTEX> &vertices){
	for(int i=0;i<vertexCount;i++){
		DeleteVertex(V);
	}
	Q_ASSERT(vertexCount==0);

	for(int i=0;i<vertices.size();i++){
		VERTEX * v=new VERTEX(vertices[i]);
		AddVertex(v);
	}
	Triangulate();
}

template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::AddVertex(VERTEX *v){
	if(vertexCount==0){
		V=v;
		lastV=v;
	}
	else{
		lastV->next=v;
		v->prev=lastV;
		lastV=v;
	}
	vertexCount++;
	VertexAdded(v);
}

template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::DeleteVertex(VERTEX *v){
	for(int i=0;i<v->Data().size();i++){
		Q_ASSERT(!v->Data()[i]->Contains((MeshVertex*)v));
	}
	if(v->prev!=nullptr){
		v->prev->next=v->next;
	}
	else{
		V=(VERTEX*)v->next;
	}
	if(v->next!=nullptr){
		v->next->prev=v->prev;
	}
	else{
		lastV=(VERTEX*)v->prev;
	}
	delete v;
	vertexCount--;
}

template<class VERTEX,class TRIANGLE>
TRIANGLE *_DynamicMesh<VERTEX, TRIANGLE>::AddTriangle(VERTEX *a, VERTEX *b, VERTEX *c){
	TRIANGLE * t=new TRIANGLE(a,b,c);
	if(triangleCount==0){
		T=t;
		lastT=t;
	}
	else{
		lastT->next=t;
		t->prev=lastT;
		lastT=t;
	}
	triangleCount++;

	a->Data().append(t);
	b->Data().append(t);
	c->Data().append(t);
	EdgeAdded(a,b);
	EdgeAdded(b,c);
	EdgeAdded(c,a);
	TriangleAdded(t);
	return t;
}

template<class VERTEX, class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::DeleteTriangle(TRIANGLE *t, bool deleteDisconnectedVertices){
	Q_ASSERT(t!=nullptr);
	for(int i=0;i<3;i++){
		VERTEX * v=(VERTEX*)(*t)(i);
		Q_ASSERT(v!=nullptr);
		v->RemoveFromTriangle((MeshTriangle*)t);
		if(deleteDisconnectedVertices && v->Data().size()==0){//vertex is not used by any triangle
			DeleteVertex(v);
		}
	}
	if(t->prev!=nullptr){
		t->prev->next=t->next;
	}
	else{
		T=(TRIANGLE*)t->next;
	}
	if(t->next!=nullptr){
		t->next->prev=t->prev;
	}
	else{
		lastT=(TRIANGLE*)t->prev;
	}
	delete t;
	triangleCount--;
}


template<class VERTEX,class TRIANGLE>
void _DynamicMesh<VERTEX,TRIANGLE>::Triangulate(){
	std::vector< Delaunay::Point > v;
	QVector<VERTEX*> vPointers;
	while(T!=nullptr){DeleteTriangle(T,false);}
	VERTEX * itVert=V;
	while(itVert!=nullptr){
		Delaunay::Point tempP;
		tempP[0] = itVert->x();
		tempP[1] = itVert->y();
		v.push_back(tempP);
		vPointers.append(itVert);
		itVert->Data().clear();
		itVert=(VERTEX*)itVert->next;
	}
	Delaunay delobject(v);
	delobject.Triangulate();
	for (Delaunay::fIterator fit = delobject.fbegin(); fit != delobject.fend(); ++fit) {
		AddTriangle(vPointers[delobject.Dest(fit)],vPointers[delobject.Org(fit)], vPointers[delobject.Apex(fit)]);
	}
}

/**** Start example implementation file ****
#include "dynamicmesh.h"
#ifndef DYNAMICMESH_IMPL_H
#define DYNAMICMESH_IMPL_H
*/
#ifndef DYNAMICMESH_IMPL_H
namespace DMesh{

class Vertex: public MeshVertex{
	public:
		Vertex():MeshVertex(){}
		Vertex(QVector2D v): MeshVertex(v){}
};

class Triangle: public MeshTriangle{
	public:
		Triangle():MeshTriangle(nullptr, nullptr, nullptr){}
		Triangle(Vertex* a, Vertex* b, Vertex* c):MeshTriangle(a,b,c){}
};

class DynamicMesh: public _DynamicMesh<Vertex,Triangle>{
	public:
		DynamicMesh():_DynamicMesh(){}
		DynamicMesh(QVector<QVector2D> p):_DynamicMesh(p){}
};

typedef _TriangleIterator<Vertex,Triangle> TriangleIterator;
typedef _VertexIterator<Vertex,Triangle> VertexIterator;
}

#endif
/**** end implementation file ****/

#endif // DYNAMICMESH_H
