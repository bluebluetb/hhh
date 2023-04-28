#include "structures.h"

bool ProblemSet::SetNextLocks(bool loop, bool guarded, double hardMinHeightFactor){ //for time series data
	//linear interpolation
	for (int i=0;i<moves.size();i++){
		QVector2D pos=moves[i].start+(moves[i].end-moves[i].start)*(1.0f/(totalIterations-1)*moveIteration);
		UnlockPoint(moves[i].v);
		if(!guarded){
			moves[i].v->setX(pos.x());
			moves[i].v->setY(pos.y());
		}
		else{
			moves[i].v->MoveTo(pos,hardMinHeightFactor);
		}
		SetLock(-moves[i].v->ID()-2,moves[i].v);
	}
	if(moveIteration==totalIterations-1 && !loop){return false;}
	moveIteration=(moveIteration+1)%totalIterations;
	return true;
}

void ProblemSet::LockMutex(QString/* debugmsg*/){
	if(!lock->tryLock()){
		quit=true;
		lock->lock();
		quit=false;
	}
	//DEBUG << "LOCK: " << debugmsg;
}
void ProblemSet::UnlockMutex(QString /*debugmsg*/){
	lock->unlock();
	//DEBUG << "UNLOCK: "<< debugmsg;
}


bool ProblemSet::Unlock(int touchId){ //returns false if it was not locked
	for(int i=0;i<TouchLocks.size();i++){
		if(TouchLocks[i].first==touchId){
			SetExtraRealWeight(TouchLocks[i].second.first,0);
			TouchLocks.removeAt(i);
			return true;
		}
	}
	return false;
}

bool ProblemSet::UnlockPoint(Vertex* vertex){ //returns false if it was not locked
	for(int i=0;i<TouchLocks.size();i++){
		if(TouchLocks[i].second.first==vertex){
			SetExtraRealWeight(TouchLocks[i].second.first,0);
			TouchLocks.removeAt(i);
			return true;
		}
	}
	return false;
}

void ProblemSet::SetExtraRealWeight(Vertex *v, double weight){
	if(!v->IsReal()){return;}
	//find real neighbours
	for(auto it=mesh->VBegin();it!=mesh->VEnd();it.Next()){//linear time is a bit slow, but you don't switch selection that often
		Vertex * v2=it.VertexPointer();
		if(v2->IsReal()){
			if(mesh->Dis(v->ID(),v2->ID())!=0){//there is actually a dissimilarity defined between the two vertices
				mesh->SetExtraW(v->ID(),v2->ID(),weight);
				if(weight>0) qDebug()<<"Setting extra selection weight for "<<v->ID()<< " - " <<v2->ID();
			}
		}
	}
}

void ProblemSet::SetLock(int touchId, Vertex *vertex){
	//is there already a lock with this touchId?
	SetExtraRealWeight(vertex,selectionWeight);
	for(int i=0;i<TouchLocks.size();i++){
		if(TouchLocks[i].first==touchId){
			TouchLocks[i].second=QPair<Vertex*,QVector2D>(vertex,(*vertex));
			return;
		}
	}
	TouchLocks.append(qMakePair(touchId,QPair<Vertex*,QVector2D>(vertex,(*vertex))));
}

Vertex* ProblemSet::HasLock(int touchId){
	for(int i=0;i<TouchLocks.size();i++){
		if(TouchLocks[i].first==touchId){
			return TouchLocks[i].second.first;
		}
	}
	return nullptr;
}

int ProblemSet::HasLockByPoint(Vertex * vertex){
	for(int i=0;i<TouchLocks.size();i++){
		if(TouchLocks[i].second.first==vertex){
			return TouchLocks[i].first;
		}
	}
	return NOLOCK;
}

void DMesh::Triangle::AddToQVector(QVector<QVector2D> &r) const{
	r.append(*(*this)(0));
	r.append(*(*this)(1));
	r.append(*(*this)(2));
}
TriangleP DMesh::Triangle::toTriangleP() const{
	QVector<QVector2D> tmp;
	AddToQVector(tmp);
	return TriangleP(tmp[0],tmp[1],tmp[2]);
}
