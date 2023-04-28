#ifndef SPARSEMAJORIZATION_H
#define SPARSEMAJORIZATION_H

#include "YaleSparseMatrix.h"
#include <QList>
#include <QPair>
#include <QVector>
#include <cmath>
#include <QDebug>
#include <iostream>
#include <functional>
#include "ProjectionMethods.h"
#include "structures.h"

#define DEBUG QDebug(QtDebugMsg)
//#define DEBUG std::cout

namespace Solver{
	void FlipTriangles(ProblemSet &p, AlgoParams &algoParams, bool guarded=true);
	void SubdivideTriangles(ProblemSet &p);
	void SubdivAndMerge(AlgoParams &params, ProblemSet &p);
}

namespace TriangleProjection
{
class SparseMajorization
{
	private:
		QVector2DContainer Z;
		QVector2DContainer tmp;
		QVector2DContainer g;
		QVector2DContainer xhat;
		int n;

		SparseMajorization(/*int m=2,*/int n):n(n){Z.resize(n);tmp.resize(n);g.resize(n);xhat.resize(n);}
		void ApplyLocksToGradient(const Vector2DContainer &x, Vector2DContainer &g);
		double innerProduct(const ArrayMapper &a, const ArrayMapper &b, int offset = 0);
		void ArrayAdd(const ArrayMapper &x, const ArrayMapper &b, ArrayMapper &r, int offset);
		double MaxAbsVal(const Vector2DContainer &x);
		void TakeDescentStep(const Vector2DContainer &x, const Vector2DContainer &d, double stepSize, Vector2DContainer &result);
		void CopyAndClear(Vector2DContainer &x, Vector2DContainer &xhat);
		double Stress(const Vector2DContainer &X) const;
		bool guardedMovement=false;
	public:
		const YaleSparseMatrix *D;
		YaleSparseMatrix *A;
		YaleSparseMatrix *B;
		int MaxIterations=10;
		bool FixFirstCoordAtZero=false;
		double lastStress;
		const double threshold = 1e-4;
		QList<QPair<int, double> > XLocks;
		QList<QPair<int, double> > YLocks;
		double maxMove;
		virtual double ProjectMethod(){return 0;}

		inline int N() const{return n;}
		static double Distance(const Vector2DContainer &x, int i, int j);

		std::function<double(int,int)> GetExtraWeight;
		double realWeight;
		SparseMajorization(BaseTask * T,const VertexContainer &X, std::function<double(int,int)> weightF=[](int, int){return 0.0;}): GetExtraWeight(weightF)
		{
			n=X.size();
			realWeight=T->algoParams.realWeight;
			Z.resize(n);tmp.resize(n);g.resize(n);xhat.resize(n);

			//D = _D + Identity(n);
			double **_D = new double*[X.size()];
			for(int i=0;i<X.size();i++){
				_D[i]=new double[X.size()];
				for(int j=0;j<X.size();j++){
					_D[i][j]=T->p.mesh->Dis(X.ID(i),X.ID(j));
				}
			}
			for(int i=0;i<n;i++){
				_D[i][i]+=1;
			}

			D = new YaleSparseMatrix(n,n,_D);
			for(int i=0;i<X.size();i++){
				delete _D[i];
			}
			delete _D;
			GetHessian(X);
		}
		virtual ~SparseMajorization(){
			delete D;
			delete A;
			delete B;
		}

		void GetHessian(const VertexContainer &X);
		//double FQ(Matrix2xN &b, Matrix2xN &x);
		void ApplyLocks(Vector2DContainer &x);
		void Gradient(Vector2DContainer &x, Vector2DContainer &b, Vector2DContainer &result, int offset = 0);

		/// <summary>
		/// optimal stepsize in direction d given gradient g:
		///  alpha = d'.g / (d'.A.d)
		/// </summary>
		/// <param name="g">gradient</param>
		/// <param name="d">descent vector</param>
		double OptimalStepSize(Vector2DContainer &g, Vector2DContainer &d, int offset = 0);
		void ReduceStress(Vector2DContainer &x);
		void Run(BaseTask*T, Vector2DContainer &x, bool flip=false);
		void GradientDescentSolve(Vector2DContainer &x, Vector2DContainer &b);
};

class ConstrainedMajorization : public SparseMajorization{
	private:
		ConstrainedProcrustesProject project;
	public:
		ConstrainedMajorization(BaseTask *T, const VertexContainer & points, ConstrainedProcrustesProject proj) : SparseMajorization(T, points,[T](int a, int b){return T->p.mesh->ExtraW(a,b);}){project=proj;}
		double ProjectMethod(){return project.Project();}
};

class WeightAdjustMajorization: public SparseMajorization{
	private:
		WeightAdjustmentProject project;
	public:
		WeightAdjustMajorization(BaseTask *T, const VertexContainer & points, WeightAdjustmentProject proj) : SparseMajorization(T, points, [T](int a, int b){return T->p.mesh->ExtraW(a,b);}){project=proj;}
		double ProjectMethod(){return project.Project();}
};
}

#endif // SPARSEMAJORIZATION_H
