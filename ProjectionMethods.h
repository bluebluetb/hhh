#ifndef PROJECTIONMETHODS_H
#define PROJECTIONMETHODS_H

#include <QVector>
#include <QVector2D>
#include <eigen3/Eigen/Core> //for Matrix
#include <eigen3/Eigen/SVD> //for SVD
#include "structures.h"
#include <climits>

typedef Eigen::Matrix2d Matrix2x2;

namespace TriangleProjection{
using namespace Eigen;

class ProjectionHelpers{
	public:
		static void ProjectXOnY(const QVector<QVector2D> &x, const QVector<QVector2D> &y, QVector<QVector2D> &result)
		{
			Q_ASSERT(x.size() == y.size());
			int n = x.size();
			QVector<QVector2D> Y;
			QVector2D m(0,0);
			for(int i=0;i<y.size();i++)
			{
				m += y[i];
			}
			m /= n;
			for (int i = 0; i < n; ++i)
			{
				Y.append(QVector2D(y[i].x() - m.x(),y[i].y() - m.y()));
			}

			//var C = X' Y;
			Matrix2x2 C;
			C(0,0)=0;
			C(1,0)=0;
			C(0,1)=0;
			C(1,1)=0;

			for (int i = 0; i < n; ++i)
			{
				double xX = x[i].x();
				double xY = x[i].y();
				double yX = Y[i].x();
				double yY = Y[i].y();
				C(0,0) += xX * yX;
				C(0,1) += xX * yY;
				C(1,0) += xY * yX;
				C(1,1) += xY * yY;
			}
			//var svd = C.Svd(true);
			Eigen::JacobiSVD<Matrix2x2> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);

			//var R = svd.VT() .TransposeThisAndMultiply(svd.U());
			Matrix2x2 VT = svd.matrixV();
			double tmp=VT(0,1);
			VT(0,1)=VT(1,0);VT(1,0)=tmp;

			Matrix2x2 U = svd.matrixU();
			Matrix2x2 R = VT * U;

			// restrict reflection (constrained procrustes)
			double det = R(0,0) * R(1,1) - R(1,0) * R(0,1);
			//var S = new DiagonalMatrix(2, 2, new double[] { 1, det });
			Matrix2x2 S = (Eigen::Matrix2d() << 1,0,0,det).finished();
			R = (VT*S) * U;

			QVector2D t;
			for (int i = 0; i < n; ++i)
			{
				//var s = xi - R.TransposeAndMultiply(yi).Transpose();
				QVector2D s = x[i] - QVector2D(R(0,0)*Y[i].x() + R(0,1)*Y[i].x(),R(1,0)*Y[i].y() + R(1,1)*Y[i].x());
				t+=s;
			}
			t /= n;
			//var TY = R.TransposeThisAndMultiply(Y.Transpose()).Transpose();
			for (int i = 0; i < n; ++i)
			{
				//var ty = TY.SubMatrix(i, 1, 0, 2);
				//result.Add(new Point(ty[0, 0] + t[0, 0], ty[0, 1] + t[0, 1]));
				result.append(QVector2D(t.x() + Y[i].x()*R(0,0)+Y[i].y()*R(0,1),
																t.y() + Y[i].x()*R(1,0)+Y[i].y()*R(1,1)));
			}
		}

		/*static QVector2D GetVector(Matrix2xN &X, int i)
		{
			return QVector2D(X(0,i), X(1,i));
		}*/

		static inline QVector2D Orthogonal(QVector2D v)
		{
			return QVector2D(-v.y(), v.x());
		}

		static QVector2D ProjectPointLine ( QVector2D a , QVector2D b , QVector2D p ) {
			return a + (b - a ) *QVector2D::dotProduct(( p - a ), (b - a )) /QVector2D::dotProduct(( b - a ),(b - a ));
		}

		static void IsFeasible(TriangleIterator &it, double factor, QVector<QVector2D> &target, bool singleMinHeight)
		{
			TriangleP t=it.GetTriangleP();
			/*bool clockwise = t.IsClockwise();

			// find smallest height of triangle
			double minHeight = std::numeric_limits<double>::infinity();
			int minimalI=-1;
			for (int i = 0; i < 3; ++i)
			{
				int j = (i + 1) % 3;
				int k = (j + 1) % 3;
				QVector2D u = t(i);
				QVector2D v = t(j);
				QVector2D w = t(k);
				QVector2D d = ProjectPointLine(u, v, w);
				QVector2D wd = w - d;
				double height = wd.length();

				if ((height < h || !clockwise) && height < minHeight)
				{
					minHeight = height;
					QVector2D o;
					o= Orthogonal(v - u);
					o.normalize();
					QVector2D neww = d - h * o;
					target.resize(3);
					target[i] = u;
					target[j] = v;
					target[k] = neww;
					minimalI=k;
				}
			}//*/
			//check if one of the heights violates the minheight of that vertex
			double minadjust=0;
			double minHeight = std::numeric_limits<double>::infinity();
			for(int i=0;i<3;i++){
				int j=(i+1)%3;
				int k=(j+1)%3;
				QVector2D u = t(i);
				QVector2D v = t(j);
				QVector2D w = t(k);
				QVector2D d = ProjectPointLine(u, v, w);
				double height = (w-d).length();
				double h;
				if(singleMinHeight) h=it.Triangle().OrigMinHeight()*factor;
				else h=it.Triangle().OrigMinHeight(i)*factor;
				if(height < h && ((!singleMinHeight && (h-height)>minadjust)||(singleMinHeight && height<minHeight))){
					QVector2D o;
					o= Orthogonal(v - u);
					o.normalize();
					QVector2D neww = d - h * o;
					target.resize(3);
					target[i] = u;
					target[j] = v;
					target[k] = neww;
					minadjust=(h-height);
					minHeight=height;
				}
			}/*
			if(minimalI==-1){
				it.Triangle().stressEdge=minimalI;
			}else{
				it.Triangle().stressEdge=(minimalI+1)%3;
			}*/
		}

		static double Project(const DMesh::Triangle &t, const QVector<QVector2D> &target, bool guarded, double minDist, bool completeOnly)
		{

			QVector<QVector2D> x;
			t.AddToQVector(x);
			QVector<QVector2D> orig=x;

			QVector<QVector2D> r;
			ProjectXOnY(x, target, r);
			double delta = 0;
			bool canProject=true;
			if(guarded && completeOnly){
				for (int i = 0; i < 3; ++i)
				{
					if(!((Vertex*)t(i))->MoveTo(r[i].x(),r[i].y(),minDist,true)){
						canProject=false;
					}
				}
			}
			if(!guarded || canProject){
				for (int i = 0; i < 3; ++i)
				{
					if(guarded){
						((Vertex*)t(i))->MoveTo(r[i].x(),r[i].y(),minDist);
					}
					else{
						t(i)->setX(r[i].x());
						t(i)->setY(r[i].y());
					}
					double d= (x[i] - orig[i]).length();
					delta +=d;
				}
			}
			else{
				qDebug() << "No complete projection possible";
			}

			return delta;
		}

		//some tests
		//[Description("X and Y match, apart from rotation, assert that projection of X on Y correctly returns a matrix matching X")]
		static void TestProcrustesProjection()
		{
			QVector<QVector2D> X;
			X.append(QVector2D(0,0));X.append(QVector2D(0,1));X.append(QVector2D(1,0));
			QVector<QVector2D> Y;
			Y.append(QVector2D(0,0));Y.append(QVector2D(1,0));Y.append(QVector2D(0,-1));

			QVector<QVector2D> P;
			ProjectionHelpers::ProjectXOnY(X, Y, P);
			//		X.Zip(P, (x, p) => AssertAreClose(x, p));
			for(int i=0;i<X.size();i++){
				Q_ASSERT((X[i]-P[i]).length()<1e-4);
			}
		}
		//[Description("Assert that projection does not allow reflection")]
		static void TestConstrainedProcrustesProjection()
		{
			QVector<QVector2D> X;
			X.append(QVector2D(0,1));X.append(QVector2D(2,0));X.append(QVector2D(3,1));
			QVector<QVector2D> Y;
			Y.append(QVector2D(0,0));Y.append(QVector2D(2,1));Y.append(QVector2D(3,0));
			QVector<QVector2D> P;
			ProjectionHelpers::ProjectXOnY(X, Y, P);
		}
};

class DummyProject{
	public:
		double Project(){return 0;}
};

class ConstrainedProcrustesProject{
	public:
		ConstrainedProcrustesProject(){}
		ConstrainedProcrustesProject(ProblemSet *p, AlgoParams &params): p(p), singleMinHeight(params.singleMinHeight),minHeightFactor(params.minProjectionHeight),hardSoftMinHeightFactor(params.softHardFactor),completeOnly(params.onlyCompleteProjections), guarded(params.guardedMovement){}
		double Project(){return ProjectTriangles();}
	private:
		ProblemSet *p=nullptr;
		bool singleMinHeight=true;
		double minHeightFactor=0.5;
		double hardSoftMinHeightFactor=0.5;
		bool completeOnly=true;
		bool guarded=true;
		double ProjectTriangles()
		{
			double delta = 0;
			for(auto it=p->mesh->VBegin();it!=p->mesh->VEnd();it.Next()){
				Vertex &v = it.Vertex();
				v.desiredPos=QVector2D(v.x(),v.y());
			}
			for(TriangleIterator it=p->mesh->Begin();it!=p->mesh->End();it.Next())
			{
				QVector<QVector2D> target;

				ProjectionHelpers::IsFeasible(it, minHeightFactor, target, singleMinHeight);
				if (target.size()>0)
				{
					delta += ProjectionHelpers::Project(it.GetTriangle(), target, guarded, minHeightFactor*hardSoftMinHeightFactor, completeOnly);
					it.Triangle().color=TGREEN;
				}else{
					it.Triangle().color=TNOCOLOR;
				}
			}
			return delta;
		}
};

class WeightAdjustmentProject{
	public:
		WeightAdjustmentProject(){}
		WeightAdjustmentProject(ProblemSet *p, AlgoParams &params): p(p), singleMinHeight(params.singleMinHeight),minHeightFactor(params.minProjectionHeight),hardSoftMinHeightFactor(params.softHardFactor),completeOnly(params.onlyCompleteProjections), guarded(params.guardedMovement){}
		double Project(){return ProjectTriangles();}
	private:
		ProblemSet *p=nullptr;
		bool singleMinHeight=true;
		double minHeightFactor=0.5;
		double hardSoftMinHeightFactor=0.5;
		bool completeOnly=true;
		bool guarded=true;
		double ProjectTriangles()
		{
			qDebug() <<"\"Projecting\" triangles";
			double delta = 0;
			for(auto it=p->mesh->VBegin();it!=p->mesh->VEnd();it.Next()){
				Vertex &v = it.Vertex();
				v.desiredPos=QVector2D(v.x(),v.y());
			}
			for(TriangleIterator it=p->mesh->Begin();it!=p->mesh->End();it.Next()){
				for(int i=0; i<3;i++){
					p->mesh->SetExtraDis(it(i).ID(),it((i+1)%3).ID(),0);
					p->mesh->SetExtraW(it(i).ID(),it((i+1)%3).ID(),0);
				}
			}
			for(TriangleIterator it=p->mesh->Begin();it!=p->mesh->End();it.Next())
			{
				double minHeight = std::numeric_limits<double>::infinity();
				int violatingPoint=-1;
				for (int i = 0; i < 3; ++i){
					int j = (i + 1) % 3;
					int k = (j + 1) % 3;
					QVector2D u = it(i);
					QVector2D v = it(j);
					QVector2D w = it(k);
					QVector2D d = ProjectPointLine(u, v, w);
					QVector2D wd = w - d;
					double height = wd.length();
					if (height < minHeight){
						minHeight = height;
						violatingPoint=k;
					}
				}
				if(violatingPoint==-1) continue;
/*
	Points / dissimilarities:
			 b3  ^
			 /|\ | h
			 \|/ v
				b
	disa/ | \ disb
		 /  |  \
		a--b2---c
	 /| disc  |\
	a2		     c2
*/
				double h=it.Triangle().OrigMinHeight()*minHeightFactor-minHeight;
				if(h<=0){continue;}
				const Vertex &a=it((violatingPoint+2)%3);
				const Vertex &b=it(violatingPoint);
				const Vertex &c=it((violatingPoint+1)%3);

				QVector2D b2=ProjectPointLine(a,c,b);
				QVector2D b3=(b-b2).normalized()*h+b;
				double b2bcAngle=atan((b2-c).length()/(b-b2).length());
				double abb2Angle=atan((a-b2).length()/(b-b2).length());
				double ddisa=(tan(abb2Angle)*(b3-b).length());
				double ddisb=(tan(b2bcAngle)*(b3-b).length());
				QVector2D a2=(a-b).normalized()*ddisa/2+a;
				QVector2D c2=(c-b).normalized()*ddisb/2+c;
				double ddisc1=(ProjectPointLine(a,c,a2)-a).length();
				double ddisc2=(ProjectPointLine(a,c,c2)-c).length();
				double ddisc=ddisc1+ddisc2;

				Q_ASSERT(it(violatingPoint).ID()==b.ID());
				p->mesh->SetExtraDis(a.ID(),b.ID(),ddisa);
				p->mesh->SetExtraDis(b.ID(),c.ID(),ddisb);
				p->mesh->SetExtraDis(c.ID(),a.ID(),-ddisc/2);
				qDebug() << "a:"<<QVector2D(a)<<"    b:"<<QVector2D(b)<<"    c:"<<QVector2D(c);
				qDebug() << "ddisa="<<ddisa<<"   ddisb="<<ddisb<<"    dddisc="<<-ddisc;
				//p->mesh->SetExtraW(a.ID(),b.ID(),100);
				//p->mesh->SetExtraW(b.ID(),c.ID(),100);
				//p->mesh->SetExtraW(c.ID(),a.ID(),1);
			}
			return delta;
		}
};
}//namespace

#endif // PROJECTIONMETHODS_H
