#define NOMINMAX

#include <limits>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "solver.h"
#include "worker.h"
#include "random.h"
#include "SparseMajorization.h"
#include "ProjectionMethods.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define MOVEMAJORS 5

using namespace TriangleProjection;

void AddLocks(SparseMajorization &maj, VertexContainer &points, ProblemSet &p){
	for(int i=0;i<p.TouchLocks.size();i++){
		for(int j=0;j<points.size();j++){
			if(points(0,j)==(*p.TouchLocks[i].second.first).x() && points(1,j)==(*p.TouchLocks[i].second.first).y()){
				maj.XLocks.append(QPair<int,double>(j,p.TouchLocks[i].second.second.x()));
				maj.YLocks.append(QPair<int,double>(j,p.TouchLocks[i].second.second.y()));
			}
		}
	}
}

void GetPoints(BaseTask * T, VertexContainer & points){
	Vertex * v=T->p.mesh->GetVertexPointer();
	for(;v!=nullptr;v=(Vertex*)v->next){
		points.append(v);
	}
}

void UpdatePoints(BaseTask *T, Vector2DContainer &points){
	Vertex * v=T->p.mesh->GetVertexPointer();
	int i=0;
	for(;v!=nullptr;v=(Vertex*)v->next){
		v->setX(points(0,i));
		v->setY(points(1,i));
		i++;
	}
	Q_ASSERT(i==points.size());
}

int Majorization(BaseTask * T) {
	T->p.LockMutex();
	VertexContainer points(T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	GetPoints(T,points);
	SparseMajorization majorization(T,points);
	AddLocks(majorization, points, T->p);
	majorization.Run(T,points);
	//UpdatePoints(T,points);
	T->info.lastStress=majorization.lastStress;
	T->p.UnlockMutex();
	return 0;
}

int _CMajorization(BaseTask * T, bool lock=true){
	if(lock){T->p.LockMutex();}
	VertexContainer points(T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	GetPoints(T,points);
	ConstrainedProcrustesProject proj(&(T->p),T->algoParams);
	ConstrainedMajorization majorization(T,points, proj);
	AddLocks(majorization, points, T->p);
	majorization.Run(T,points);
	//UpdatePoints(T,points);
	T->info.lastStress=majorization.lastStress;
	if(lock){T->p.UnlockMutex();}
	return 0;
}

int CMajorization(BaseTask* T){
	return _CMajorization(T);
}

int WMajorization(BaseTask * T){
	T->p.LockMutex();
	VertexContainer points(T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	GetPoints(T,points);
	WeightAdjustmentProject proj(&(T->p),T->algoParams);
	int limit=T->algoParams.limit;
	T->algoParams.limit=1;
	WeightAdjustMajorization majorization(T,points, proj);
	AddLocks(majorization, points, T->p);
	//backup values
	for(auto it=T->p.mesh->VBegin();it!=T->p.mesh->VEnd();it.Next()){
		it.VertexPointer()->Backup();
	}
	majorization.Run(T,points); //adjusts the weights
	//revert
	for(auto it=T->p.mesh->VBegin();it!=T->p.mesh->VEnd();it.Next()){
		it.VertexPointer()->Revert();
	}
	SparseMajorization major2(T,points);
	AddLocks(major2,points,T->p);
	major2.Run(T,points);	//now we do a normal majorization with the adjusted weights
	T->algoParams.limit=limit;
	T->info.lastStress=majorization.lastStress;
	T->p.UnlockMutex();
	return 0;
}

int CSubDivMajor(BaseTask* T){
	for(int i=0;i<10;i++){_CMajorization(T);}
	//Solver::SubdivideTriangles(T->p);
	Solver::SubdivAndMerge(T->algoParams,T->p);
	return 0;
}

int CFixSubDivMajor(BaseTask* T){
	_CMajorization(T);
	Solver::FlipTriangles(T->p,T->algoParams);
//	Solver::SubdivideTriangles(T->p);
	Solver::SubdivAndMerge(T->algoParams,T->p);
	return 0;
}

int CSubDivFixMajor(BaseTask* T){
	_CMajorization(T);
	//Solver::SubdivideTriangles(T->p);
	Solver::SubdivAndMerge(T->algoParams,T->p);
	Solver::FlipTriangles(T->p,T->algoParams);
	return 0;
}

int MoveCMajoriazation(BaseTask * T){
	T->p.LockMutex();
	T->p.SetNextLocks(false,T->algoParams.guardedMovement,T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	for(int i=0;i<MOVEMAJORS;i++){
		_CMajorization(T,false);
	}
	T->p.UnlockMutex();
	return 0;
}

int MoveSubdivCMajor(BaseTask * T){
	T->p.LockMutex();
	T->p.SetNextLocks(false,T->algoParams.guardedMovement,T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	for(int i=0;i<MOVEMAJORS;i++){
		_CMajorization(T,false);
		//Solver::SubdivideTriangles(T->p);
		Solver::SubdivAndMerge(T->algoParams,T->p);
	}
	T->p.UnlockMutex();
	return 0;
}

int MoveCMajoriazationL(BaseTask * T){
	T->p.LockMutex();
	T->p.SetNextLocks(true,T->algoParams.guardedMovement,T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	for(int i=0;i<MOVEMAJORS;i++){
		_CMajorization(T,false);
	}
	T->p.UnlockMutex();
	return 0;
}

int Test(BaseTask *){
	ProjectionHelpers::TestProcrustesProjection();
	ProjectionHelpers::TestConstrainedProcrustesProjection();
	return 0;
}

int FixMajorization(BaseTask *T){
	T->p.LockMutex();
	_CMajorization(T,false);
	Solver::FlipTriangles(T->p,T->algoParams);
	T->p.UnlockMutex();
	return 0;
}

int FixMoveMajor(BaseTask* T){
	T->p.LockMutex();
	T->p.totalIterations=100;
	T->p.SetNextLocks(false,T->algoParams.guardedMovement,T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	for(int i=0;i<5;i++){
		_CMajorization(T,false);
		Solver::FlipTriangles(T->p,T->algoParams);
	}
	T->p.UnlockMutex();
	T->p.UnlockMutex();
	return 0;
}

int FixMoveMajor2(BaseTask* T){
	T->p.LockMutex();
	T->p.totalIterations=100;
	T->p.SetNextLocks(false,T->algoParams.guardedMovement,T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	VertexContainer points(T->algoParams.minProjectionHeight*T->algoParams.softHardFactor);
	GetPoints(T,points);
	for(int i=0;i<5;i++){
		ConstrainedProcrustesProject proj(&(T->p),T->algoParams);
		ConstrainedMajorization majorization(T, points, proj);
		AddLocks(majorization, points, T->p);
		majorization.Run(T,points,true);
		//UpdatePoints(T,points);
		T->info.lastStress=majorization.lastStress;
	}
	T->p.UnlockMutex();
	T->p.UnlockMutex();
	return 0;
}

int Fixer(BaseTask *T){
	T->p.LockMutex();
	Solver::FlipTriangles(T->p,T->algoParams);
	T->p.UnlockMutex();
	return 0;
}

Solver::type Solver::solvers[] = {
	{"Majorization", Majorization }, //0
	{"CMajorization", CMajorization }, //1
	{"WMajorization", WMajorization }, //2
	{"MoveCMajor", MoveCMajoriazation}, //3
	{"MoveCMajorL", MoveCMajoriazationL}, //4
	{"FixCMajor", FixMajorization},//5
	{"MoveFixCMajor", FixMoveMajor},//6
	{"MoveFixCMajor2", FixMoveMajor2},//7
	{"SubdivCMajor", CSubDivMajor}, //8
	{"MoveSubCMajor",MoveSubdivCMajor},//9
	{"SubdivFixCMajor",CSubDivFixMajor},//10
	{"FixSubdivCMajor",CFixSubDivMajor},//11
	{"Fixer",Fixer},//12
	//{"Test", Test }, //13
	{nullptr,nullptr}, // Last element of list.
};
