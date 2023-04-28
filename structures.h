#ifndef STRUCTURES_H
#define STRUCTURES_H

#ifndef nullptr
#define nullptr NULL
#endif

#include<QVector>
#include<QVector2D>
#include<QVector3D>
#include<QDebug>
#include<QMutex>
#include<limits>
#include"DynamicMesh_impl.h"

#define DEBUG QDebug(QtDebugMsg)
#define NOLOCK std::numeric_limits<int>::min()

inline double min(double a, double b){
	return (a<b)?a:b;
}

inline double max(double a, double b){
	return (a>b)?a:b;
}



/*
class ArrayMapper{
	public:
		ArrayMapper(): p(nullptr),mesh(nullptr),which(0){}
		ArrayMapper(DynamicMesh *mesh, unsigned int i): mesh(mesh),p(mesh->GetVertexPointer()),which(i){}
		void Set(double val){return (which==0)?p->setX(val):p->setY(val);}
		double Val() const{return (which==0)?p->x():p->y();}

		int size()const{return mesh->nPoints();}
		void Next(){if(p!=nullptr) p=(Vertex*)p->next;}
		bool AtEnd(){return (p==nullptr);}
		void Reset(){p=mesh->GetVertexPointer();}

		double Get(int i){  //expensive linear operation from when there was still an array of vertices, trying to remove it from the code
			Vertex* backup=p;
			Reset();
			for(int j=0;j<i;j++){
				Next();
			}
			double ret=Val();
			p=backup;
			return ret;
		}
	private:
		DynamicMesh * mesh;
		Vertex *p;
		unsigned int which;
};*/

#include <GL/gl.h>
struct MovementData{
		MovementData(): v(nullptr){}
		MovementData(Vertex*  v, QVector2D start, QVector2D end):v(v), start(start), end(end){}
		Vertex * v;
		QVector2D start;
		QVector2D end;
};

class ProblemSet{
	public:
		ProblemSet(){/*XLocks=QList<QPair<int,double> >();YLocks=QList<QPair<int,double> >();*/lock=new QMutex();quit=false;mesh=new DynamicMesh();}
		//ProblemSet(Matrix2xN _X, double ** _D) : X(_X) , D(_D), lock(new QMutex()),quit(false) {}
		int size(){return mesh->nPoints();}
		//void AddLock(int index, double x, double y){XLocks.append(QPair<int,double>(index,x));YLocks.append(QPair<int,double>(index,y));}
		//QList<QPair<int,double> > XLocks;
		//QList<QPair<int,double> > YLocks;
		void Init(int size); //allocates and initializes D, W and X
		QVector2D GetErrorVec(int a);

		//<touchId (-1 for normal mouse), <pointId, Position> >. touchIds -MaxInt to -2 are used for extra locks of point -(id+2)
		QList<QPair<int, QPair<Vertex*,QVector2D> > > TouchLocks;

		QVector<MovementData> moves;
		int moveIteration=0;
		int totalIterations=100;
		bool SetNextLocks(bool loop, bool guarded=false, double hardMinHeightFactor=0.25); //for time series data, returns false iff loop==false and moveIteration==totalIterations-1
		bool Unlock(int touchId); //returns false if it was not locked
		bool UnlockPoint(Vertex* vertex); //returns false if it was not locked
		void SetLock(int touchId, Vertex* vertex);
		Vertex *HasLock(int touchId); //returns nullptr if no lock, Pointer otherwise
		int HasLockByPoint(Vertex *vertex); //returns touchId if locked, NOLOCK otherwise
		//bool isClockwiseTriangle(TriangleProjection::Triangle t);
		bool hasAssociatedValues=false; //do the real vertices have a value attribute?

		DMesh::DynamicMesh * mesh;
		GLuint glTextureId;
		bool hasTexture;
		QVector2D minOrigCoord=QVector2D(0,0);
		QVector2D maxOrigCoord=QVector2D(0,0);
		double minTriangleArea=240;
		int imgWidth=0,imgHeight=0;
		void SetExtraRealWeight(Vertex* v,double weight);
		double selectionWeight=4;
	private:
		QMutex *lock; //toggled by the solver and the glwidget
	public:
		volatile bool quit;
		void LockMutex(QString debugmsg="");
		void UnlockMutex(QString debugmsg="");
};

struct TaskInfo{
		double executionTime;
		double lastStress;
		static int solves;
		TaskInfo(): executionTime(0),lastStress(0){}
};

struct AlgoParams{
		AlgoParams(int algoID=0, int limit=0, double minProjectionHeight=0.5, double softHardFactor=0.5, double realWeight=1, bool onlyCompleteProjections=false, bool guardedMovement=true, int solvesPerSec=30, bool continuousSolve=false): algoID(algoID), limit(limit), minProjectionHeight(minProjectionHeight), softHardFactor(softHardFactor), realWeight(realWeight),  onlyCompleteProjections(onlyCompleteProjections), guardedMovement(guardedMovement), solvesPerSec(solvesPerSec),continuousSolve(continuousSolve){}
		int algoID;
		int limit; //maximum number of stress reduce executions per solve, 0=no limit
		double minProjectionHeight;
		double maxProjectionHeight=3;
		double softHardFactor;
		double realWeight;
		bool onlyCompleteProjections;
		bool guardedMovement;
		int solvesPerSec;
		bool singleMinHeight=true;
		bool continuousSolve=false;
};

class BaseTask {
	public:
		TaskInfo info;
		BaseTask(ProblemSet &p, AlgoParams &algoParams) : p(p), algoParams(algoParams){}
		virtual ~BaseTask(){}
		ProblemSet &p;
		AlgoParams &algoParams;
};

struct GenParams{
		GenParams(int genID=0, QString image="", int count=12, int maxc=100, int seed=42, int grid=10) : genID(genID), image(image),count(count),maxc(maxc),seed(seed), grid(grid){}
		int genID;
		QString image;
		int count;
		int maxc;
		int seed;
		int grid;
		bool overlayGrid=true;
		int overlayGridX=10;
		int overlayGridY=10;
};

struct VisParams{
		VisParams(bool showPoints=true, bool showRealPoints=true, bool showLockedPoints=true, bool showSelectedPoints=true, bool showEdges=true, bool showIso=false, bool smoothPoints=true, int lineWidth=1, int pointSize=6, int isoGrid=30, bool colors=false, bool autoSolve=true, double arrowWeight=0):
			showPoints(showPoints),showRealPoints(showRealPoints), showLockedPoints(showLockedPoints), showSelectedPoints(showSelectedPoints), showEdges(showEdges), showIso(showIso), smoothPoints(smoothPoints), lineWidth(lineWidth), pointSize(pointSize), isoGrid(isoGrid), colors(colors), autoSolve(autoSolve), arrowWeight(arrowWeight){}
		bool showPoints;
		bool showRealPoints;
		bool showLockedPoints;
		bool showSelectedPoints;
		bool showEdges;
		bool showIso;
		bool smoothPoints;
		int lineWidth;
		int pointSize;
		int isoGrid;
		bool colors;
		bool autoSolve;
		double arrowWeight;
		bool bridgesRidges=true;
		double bridgeFactor=0.2;//if more than this factor of, draw bridge/ridge
};


#endif // STRUCTURES_H
