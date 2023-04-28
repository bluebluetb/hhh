#include "SparseMajorization.h"

#define STEPFACTOR 1

using namespace TriangleProjection;

#include <QVector3D>
#include <cmath>

/*inline QVector2D ProjectPointLine ( QVector2D a , QVector2D b , QVector2D p ) {
	return a + (b - a ) *QVector2D::dotProduct(( p - a ), (b - a )) /QVector2D::dotProduct(( b - a ),(b - a ));
}*/

void Solver::SubdivAndMerge(AlgoParams &params, ProblemSet &p){
	//p.LockMutex();
	for(auto it=p.mesh->Begin();it!=p.mesh->End();it.Next()){it.Triangle().changed=false;}
	for(auto it=p.mesh->VBegin();it!=p.mesh->VEnd();it.Next()){
		Vertex &v=it.Vertex();
		if(!v.IsReal()){continue;}
		if(true||(v-v.desiredPos).lengthSquared()>0.01){
			//last procrustes proj was unable to place point at desired position-> can we fix this through subdiv/merging?
			//find restricting triangles
			bool processed=false;
			do{
				processed=false;
				for(int t=0;t<v.Data().size();t++){
				auto it=p.mesh->Iterator((Triangle*)v.Data()[t]);
				if(it.Triangle().changed){continue;}
				Vertex * one, * two;
				int vi=-1;
				for(int i=0;i<3;i++){
					if(it(i)==v){
						one=&it[(i+1)%3];
						two=&it[(i+2)%3];
						vi=i;break;
					}
				}if(vi==-1){
					qDebug() << "blaat";
				}
				//check if this triangle prevents movement in our desired direction
				QVector2D proj=ProjectPointLine(*one, *two, v);
				double d1=(v-proj).length();
				double mh;
				if(params.singleMinHeight) mh=it.Triangle().OrigMinHeight(vi);
				else mh=it.Triangle().OrigMinHeight(vi);
				mh*=params.minProjectionHeight;
				if(d1<mh*0.9){
					VertexIterator vIt;
					if(!one->IsReal() || !two->IsReal()){
						if(((*one)-v).lengthSquared() < ((*two)-v).lengthSquared()){
							if(!one->IsReal()){vIt=p.mesh->VIterator(one);}
							else{vIt=p.mesh->VIterator(two);}
						}
						else{
							if(!two->IsReal()){vIt=p.mesh->VIterator(two);}
							else{vIt=p.mesh->VIterator(one);}
						}
						vIt.RemoveAndMerge();
					}
					else{
						QVector2D ot=(*one-(*two));
						QVector2D vt=(v-(*two));
						if(!(0 < QVector2D::dotProduct(vt,ot) && QVector2D::dotProduct(vt,ot)< QVector2D::dotProduct(ot,ot))){
							continue;
						}
						it.SubDivide((vi+1)%3,
												 [&p,&proj](Vertex * , Vertex * , Triangle * ){
							Vertex * v=new Vertex(proj);
							v->SetOrigCoords();//for disimilarity calculations
							//QVector2D loc=(((*a)+(*b))/2);
							QVector2D loc=proj;//ProjectPointLine(*a,*b,*c);
							v->setX(loc.x());
							v->setY(loc.y());
							double w=p.maxOrigCoord.x()-p.minOrigCoord.x();
							double h=p.maxOrigCoord.y()-p.minOrigCoord.y();
							QVector2D x=v->GetOrigCoords();
							v->tx=x.x()/w;
							v->ty=x.y()/h; //for texure mapping
							return v;
						}
						);
					}
					it.Triangle().changed=true;
					processed=true;
				}
				}
			}while(processed==true);//continue as long as we are still making some progress
		}
	}
	//p.UnlockMutex();
}

/*void Solver::SubdivideTriangles(ProblemSet $p){
	for(auto it=p.mesh->Begin();it!=p.mesh->End();it.Next()){
		if(it.Triangle().stressEdge!=-1){
			TriangleP t(it[0].GetOrigCoords(),it[1].GetOrigCoords(),it[2].GetOrigCoords());
			if(t.Area()>p.minTriangleArea){
				it.SubDivide(it.Triangle().stressEdge,
					[&p](Vertex * a, Vertex * b, Triangle * t){
						//QVector2D splitPos((a->GetOrigCoords()+b->GetOrigCoords())/2); //TODO: do something smarter here
						Vertex* c;
						for(int i=0;i<3;i++){if((*t)(i)!=a && (*t)(i)!=b){c=(Vertex*)(*t)(i);}}
						QVector2D splitPos=ProjectPointLine(*a,*b,*c);//ProjectPointLine(a->GetOrigCoords(),b->GetOrigCoords(),c->GetOrigCoords());
						Vertex * v=new Vertex(splitPos);
						v->SetOrigCoords();//for disimilarity calculations
						//QVector2D loc=(((*a)+(*b))/2);
						QVector2D loc=splitPos;ProjectPointLine(*a,*b,*c);
						v->setX(loc.x());
						v->setY(loc.y());
						double w=p.maxOrigCoord.x()-p.minOrigCoord.x();
						double h=p.maxOrigCoord.y()-p.minOrigCoord.y();
						QVector2D x=v->GetOrigCoords();
						v->tx=x.x()/w;
						v->ty=x.y()/h; //for texure mapping
						return v;
					}
				);
			}
			it.Triangle().stressEdge=-1;
		}
	}
}*/

void Solver::FlipTriangles(ProblemSet &p, AlgoParams &algoParams, bool guarded){
	TriangleIterator it;
	for(it=p.mesh->Begin();it!=p.mesh->End();it.Next()){
		if(it.Triangle().color==TGREEN || it.Triangle().color==TBLUE){
			//triangle might benefit from a flip
			int longestEdge=-1;
			double len2=-1;
			for(int i=0;i<3;i++){
				double l=it.GetEdge(i).LengthSquared();
				if(l>len2){longestEdge=i;len2=l;}
			}
			if(!it.HasNeighbour(longestEdge)){return;}
			//if(sqrt(len2)>algoParams.minProjectionHeight*2){//flip would definately result in feasible triangles
				Vertex * stressV=nullptr;
				if(it.Triangle().stressEdge!=-1){stressV=&it[it.Triangle().stressEdge];}
				it.Triangle().color=TNOCOLOR;
				TriangleIterator it2=it.GetNeighbourIt(longestEdge);
				it2.Triangle().color=TNOCOLOR;
				QPair<QVector2D,QVector2D> edgeV=qMakePair(QVector2D(it[longestEdge].GetOrigCoords()),QVector2D(it[(longestEdge+1)%3].GetOrigCoords()));
				it.FlipEdge(longestEdge);
				it.Triangle().stressEdge=-1;
				TriangleP origT(it[0].GetOrigCoords(),it[1].GetOrigCoords(),it[2].GetOrigCoords());
				TriangleP origT2(it2[0].GetOrigCoords(),it2[1].GetOrigCoords(),it[3].GetOrigCoords());
				//if(guarded && ((origT.Encloses(edgeV.first) && !origT.Contains(edgeV.first)) || (origT2.Encloses(edgeV.second) && !origT2.Contains(edgeV.second)))){
				if(guarded && ((origT.Encloses(edgeV.first)) || (origT2.Encloses(edgeV.second)))){
					it.FlipEdge((longestEdge+2)%3);//revert flip
					it2.Triangle().color=TBLUE;
					//fix stressv
					for(int i=0;i<3;i++){if(&it[i]==stressV) it.Triangle().stressEdge=i;}
				}
				else{
					if(it.GetTriangleP().MinHeight()<algoParams.minProjectionHeight || it2.GetTriangleP().MinHeight()<algoParams.minProjectionHeight){
						//flip resulted in incorrect triangles, revert
						it.FlipEdge((longestEdge+2)%3);//revert flip
						//fix stressv
						for(int i=0;i<3;i++){if(&it[i]==stressV) it.Triangle().stressEdge=i;}
					}
				}
			//}
		}
	}
}

void SparseMajorization::GetHessian(const VertexContainer &X)
{
	A= new YaleSparseMatrix(D->A, D->IA, D->JA);
	B= new YaleSparseMatrix(D->A, D->IA, D->JA);

	for (int i = 0; i < n; ++i)
	{
		int endRow = A->IA[i + 1];
		int ii = 0;
		double vii = 0;
		for (int row = A->IA[i]; row < endRow; ++row)
		{
			int j = A->JA[row];
			if (i == j)
			{
				ii = row;
			}
			else
			{
				double d = A->A[row];
				double w = 1;
				if(X.IsReal(i) && X.IsReal(j)){
					w=realWeight;
				}
				w+=GetExtraWeight(X.ID(i),X.ID(j));
				w/=d;
				w *= -w;
				A->A[row] = w;
				vii -= w;
			}
		}
		A->A[ii] = vii;
	}

}
/*
double SparseMajorization::FQ(Matrix2xN &b, Matrix2xN &x)
{
	int m = b.size();
	Q_ASSERT(m == x.size());
	int n = b.rows();
	Q_ASSERT(n == x.rows());
	double xb = 0;
	double xAx = 0;
	for (int i = 0; i < m; ++i)
	{
		xb += innerProduct(x[i], b[i]);
		A->MultiplyByVector(x[i], tmp[0]);
		xAx += innerProduct(x[i], tmp[0]);
	}
	return 0.5 * xAx + xb;
}*/

void SparseMajorization::ApplyLocks(Vector2DContainer &x)
{
	QList<QPair<int, double> >::const_iterator  it=XLocks.begin();
	while(it!=XLocks.constEnd())
	{
		x[0].Set(it->first,it->second);
		++it;
	}

	it=YLocks.begin();
	while(it!=YLocks.constEnd())
	{
		x[1].Set(it->first,it->second);
		++it;
	}
}

void SparseMajorization::Gradient(Vector2DContainer &x, Vector2DContainer &b, Vector2DContainer &result, int offset)
{
	int n = x.size();
	Q_ASSERT(n == b.size());
	Q_ASSERT(n == result.size());
	for (int i = 0; i < offset; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			result[j].Set(i,0);
		}
	}
	for (int i = 0; i < 2; ++i)
	{
		A->MultiplyByVector(x[i], result[i], offset);
		ArrayAdd(result[i], b[i], result[i], offset);
	}
}

double SparseMajorization::OptimalStepSize(Vector2DContainer &g, Vector2DContainer &d, int offset)
{
	int m = 2;
	int n = g.size();
	Q_ASSERT(n == d.size());
	double numerator = 0;
	double denominator = 0;
	for (int i = 0; i < m; ++i)
	{
		numerator += innerProduct(g[i], d[i], offset);
		A->MultiplyByVector(d[i], tmp[0], offset);
		denominator += innerProduct(d[i], tmp[0], offset);
	}
	return (numerator / denominator) * STEPFACTOR;
}

//private stuff//

double SparseMajorization::Stress(const Vector2DContainer &X) const
{
	double sum = 0;
	for (int i = 0; i < n; ++i)
	{
		int endRow = A->IA[i + 1];
		for (int row = A->IA[i]; row < endRow; ++row)
		{
			int j = A->JA[row];
			if (i == j) continue;
			double dx = X(0,i) - X(0,j);
			double dy = X(1,i) - X(1,j);
			double dist = sqrt(dx * dx + dy * dy);
			double diff = dist - D->A[row];
			sum += -diff * diff * A->A[row];
		}
	}
	return sum;
}

void SparseMajorization::ApplyLocksToGradient(const Vector2DContainer &x, Vector2DContainer &g)
{
	QList<QPair<int, double> >::const_iterator it=XLocks.begin();
	while(it!=XLocks.constEnd())
	{
		double p = it->second;
		int j = it->first;
		g[0].Set(j,x(0,j) - p);
		++it;
	}

	it=YLocks.begin();
	while(it!=YLocks.constEnd())
	{
		double p = it->second;
		int j = it->first;
		g[1].Set(j,x(1,j) - p);
		++it;
	}

}

double SparseMajorization::innerProduct(const ArrayMapper &a, const ArrayMapper &b, int offset)
{
	int n = a.size();
	double sum = 0;
	for (int i = offset; i < n; ++i)
	{
		sum += a(i) * b(i);
	}
	return sum;
}


void SparseMajorization::ArrayAdd(const ArrayMapper &x, const ArrayMapper &b, ArrayMapper &r, int offset)
{
	for (int i = offset; i < x.size(); ++i)
	{
		r.Set(i,x(i) + b(i));
	}
}

double SparseMajorization::MaxAbsVal(const Vector2DContainer &x)
{
	int m = 2;
	int n = x.size();
	double max = 0;
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			double a = fabs(x(i,j));
			if (a > max) max = a;
		}
	}
	return max;
}

void SparseMajorization::TakeDescentStep(const Vector2DContainer &x, const Vector2DContainer &d, double stepSize, Vector2DContainer &result)
{
	int m = 2; //m=2
	int n = x.size();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			result[i].Set(j,x(i,j) - stepSize * d(i,j));
		}
	}
}

void SparseMajorization::CopyAndClear(Vector2DContainer &x, Vector2DContainer &xhat)
{
	int n = x.size();
	for (int i = 0; i < n; ++i)
	{
		if(!guardedMovement){
			x.Set(i,xhat(0,i),xhat(1,i));
		}else{
			x.MoveTo(i,xhat(0,i),xhat(1,i));
		}
		xhat.Set(i,0,0);
	}
}

double SparseMajorization::Distance(const Vector2DContainer &x, int i, int j)
{
	double xi = x(0,i), yi = x(1,i);
	double xj = x(0,j), yj = x(1,j);
	double dx = xi - xj;
	double dy = yi - yj;
	return sqrt(dx * dx + dy * dy);
}


void SparseMajorization::ReduceStress(Vector2DContainer &x)
{
	for (int i = 0; i < A->RowCount(); ++i)
	{
		int endRow = A->IA[i + 1];
		int ii = 0;
		double d = 0;
		for (int row = A->IA[i]; row < endRow; ++row)
		{
			int j = A->JA[row];
			if (i == j)
			{
				ii = row;
				continue;
			}
			double a_ij = A->A[row];
			double d_ij = D->A[row];
			double dist = Distance(x, i, j);
			d -= B->A[row] = -a_ij * d_ij / dist;
		}
		B->A[ii] = d;
	}
	B->MultiplyByVector(x[0], Z[0]);
	B->MultiplyByVector(x[1], Z[1]);
	GradientDescentSolve(x, Z);
}

void SparseMajorization::Run(BaseTask *T, Vector2DContainer &X, bool flip)
{
	int limit=T->algoParams.limit;
	guardedMovement=T->algoParams.guardedMovement;
	Q_ASSERT(X.size() == n);
	lastStress = Stress(X);
	//DEBUG << "init stress = " << lastStress;
	double reductionRatio;
	int runs=0;
	do
	{
		qDebug()<<"Running \"ReduceStress(X)\"";
		ReduceStress(X);
		double stress = Stress(X);
		reductionRatio = (lastStress - stress) / lastStress;
		lastStress = stress;
		//DEBUG <<"stress = " << stress << "\treduction ratio = "<<reductionRatio;
		runs++;
		if(flip){Solver::FlipTriangles(T->p,T->algoParams);}
		if(T->p.quit){return;}
	} while (lastStress > threshold && reductionRatio > threshold && runs!=limit);
	//DEBUG << "LastStress= " << lastStress;
}

void SparseMajorization::GradientDescentSolve(Vector2DContainer &x, Vector2DContainer &b)
{
	//DEBUG << "SparseGradientDescentSolve, Finit = " << FQ(b, x);
	int m = x.size(); //m=2
	int offset = FixFirstCoordAtZero ? 1 : 0;
	for (int i = 0; i < MaxIterations; ++i)
	{
		if (FixFirstCoordAtZero)
		{
			for (int j = 0; j < m; ++j)
			{
				x[j].Set(0,0);
			}
		}
		ApplyLocks(x);
		Gradient(x, b, g, offset);
		ApplyLocksToGradient(x, g);
		if (MaxAbsVal(g) < 1e-3)
		{
			ProjectMethod();
			break;
		}
		double a = OptimalStepSize(g, g, offset);
		if (std::isnan(a) || a < 1e-2)
		{
			ProjectMethod();
			break;
		}
		TakeDescentStep(x, g, a, xhat);
		CopyAndClear(x, xhat);

		double delta;
		int projectionCount = 0;
		do
		{
			delta = 0;
			//ApplyLocks(x);
			delta += ProjectMethod();
		}
		while (delta > 1 && ++projectionCount < 1);
		//DEBUG << "SparseGradientDescentSolve, F["<<i<<"] = "<< FQ(b, x);
	}
}

