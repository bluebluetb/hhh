#ifndef QUADTREE_MESH_H
#define QUADTREE_MESH_H

//quick and dirty (and inefficient) approach to creating the points for a non-uniform quadtree based mesh

#include <QVector>
#include <QVector2D>
#include <QRectF>

void AddPoints(QVector<QVector2D> &meshpts, QRectF rect,int extraDivisions=0){
	if(extraDivisions>=2){
		meshpts.append(QVector2D(rect.topLeft()));
		meshpts.append(QVector2D(rect.topRight()));
		meshpts.append(QVector2D(rect.bottomLeft()));
		meshpts.append(QVector2D(rect.bottomRight()));
	}else{
		QRectF NWr(rect.topLeft(),QPointF(rect.left()+rect.width()/2,rect.top()+rect.height()/2));
		QRectF NEr(NWr.topRight(),NWr.size());
		QRectF SWr(NWr.bottomLeft(),NWr.size());
		QRectF SEr(NWr.bottomRight(),NWr.size());
		AddPoints(meshpts,NWr,extraDivisions+1);
		AddPoints(meshpts,NEr,extraDivisions+1);
		AddPoints(meshpts,SWr,extraDivisions+1);
		AddPoints(meshpts,SEr,extraDivisions+1);
	}
}

void GetPoints(QVector<QVector2D> &pts, QRectF rect, QVector<QVector2D> &meshpts, std::function<bool(QVector<QVector2D> &,QRectF &)> canBranch){
	if(canBranch(pts,rect)){
		QVector<QVector2D> NW,NE,SW,SE;
		QRectF NWr(rect.topLeft(),QPointF(rect.left()+rect.width()/2,rect.top()+rect.height()/2));
		QRectF NEr(NWr.topRight(),NWr.size());
		QRectF SWr(NWr.bottomLeft(),NWr.size());
		QRectF SEr(NWr.bottomRight(),NWr.size());
		for(int i=0;i<pts.size();i++){
			if(NWr.contains(pts[i].toPointF())){
				NW.append(pts[i]);
				continue;
			}
			if(NEr.contains(pts[i].toPointF())){
				NE.append(pts[i]);
				continue;
			}
			if(SWr.contains(pts[i].toPointF())){
				SW.append(pts[i]);
				continue;
			}
			if(SEr.contains(pts[i].toPointF())){
				SE.append(pts[i]);
				continue;
			}
		}
		GetPoints(NW,NWr,meshpts,canBranch);
		GetPoints(NE,NEr,meshpts,canBranch);
		GetPoints(SW,SWr,meshpts,canBranch);
		GetPoints(SE,SEr,meshpts,canBranch);
	}else{
		AddPoints(meshpts,rect);
	}
}


void MeshPoints(QVector<QVector2D> &realpts, QRectF rect, QVector<Vertex> &meshpts, double maxSide){
	QVector<QVector2D> tmp;
	GetPoints(realpts,rect,tmp,
		[&](QVector<QVector2D> &p,QRectF& r){
		if(r.height()>maxSide){return true;}
		if(p.size()>1){return true;}
		return false;
	});
	//remove duplicates
	meshpts.append(tmp[0]);
	for(int i=1;i<tmp.size();i++){
		bool duplicate=false;
		for(int j=0;j<meshpts.size();j++){
			if(tmp[i]==meshpts[j]){
			 duplicate=true; break;
			}
		}
		if(!duplicate){meshpts.append(tmp[i]);}
	}
}

#endif // QUADTREE_MESH_H
