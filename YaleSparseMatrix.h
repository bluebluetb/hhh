#ifndef YALESPARSEMATRIX_H
#define YALESPARSEMATRIX_H

#include <cstdio>
#include <QList>
#include <QVector>
#include "structures.h"

namespace TriangleProjection
{
/// <summary>
/// The Yale Sparse Matrix Format stores an initial sparse m×n matrix, M,
/// in row form using three one-dimensional arrays. Let NNZ denote the number
/// of nonzero entries of M. The first array is A, which is of length NNZ, and
/// holds all nonzero entries of M in left-to-right top-to-bottom (row-major) order.
/// The second array is IA, which is of length  (i.e., one entry per row, plus one).
/// IA(i) contains the index in A of the first nonzero element of row i.
/// Row i of the original matrix extends from A(IA(i)) to A(IA(i+1)-1),
/// i.e. from the start of one row to the last index before the start of the next.
/// The third array, JA, contains the column index (zero-based) of each element of A,
/// so it also is of length NNZ. For example, the matrix
///
/// [ 1 2 0 0 ]
/// [ 0 3 9 0 ]
/// [ 0 1 4 0 ]
///
/// is a three-by-four matrix with six nonzero elements, so
///
/// A = [ 1 2 3 9 1 4 ]   (array of non-zero element values)
/// IA = [ 0 2 4 6 ]   (array of index of first nonzero element of row i)
/// JA = [ 0 1 1 2 1 2 ]   (array of column index of each A element)
/// </summary>
class YaleSparseMatrix
{
	private:
		int rows;
	public:
		/// <summary>
		/// array of non-zero element values
		/// </summary>
		QVector<double> A;

		/// <summary>
		/// array of index of first nonzero element of row i
		/// </summary>
		QVector<unsigned> IA;

		/// <summary>
		/// array of column index of each A element
		/// </summary>
		QVector<unsigned> JA;

		int RowCount()
		{
			return rows;
		}

		YaleSparseMatrix(const QVector<double> &_A,const QVector<unsigned> &_IA, const QVector<unsigned> &_JA) : A(_A), IA(_IA),JA(_JA)
		{
			Q_ASSERT(A.size() == JA.size());
			rows =IA.size() -1;
		}

		YaleSparseMatrix(int _rows, int cols, double** dense): rows(_rows)
		{
			QList<double> lA;
			IA.resize(rows + 1);
			QList<unsigned int> lJA;
			for (int i = 0; i < rows; ++i)
			{
				IA[i] = lA.size();
				bool first = true;
				for (int j = 0; j < cols; ++j)
				{
					double v = dense[i][j];
					if (v != 0)
					{
						if (first)
						{
							IA[i] = lA.size();
							first = false;
						}
						lA.append(v);
						lJA.append(j);
					}
				}
			}
			IA[rows] = lA.size();
			this->A =  lA.toVector();
			this->JA = lJA.toVector();
		}

		void MultiplyByVector(const ArrayMapper &x, ArrayMapper & r, int offset = 0)
		{
			Q_ASSERT(x.size() == rows);
			Q_ASSERT(r.size() == rows);
			Q_ASSERT(rows == IA.size() - 1);
			for (int i = offset; i < rows; ++i)
			{
				int endRow = IA[i + 1];
				r.Set(i,0);
				for (int row = IA[i]; row < endRow; ++row)
				{
					int j = JA[row];
					if (j < offset) continue;
					double a_ij = A[row];
					double val=r(i)+a_ij * x(j);
					r.Set(i,val);
				}
			}
		}

		//not used, so not checked!
		YaleSparseMatrix MultiplyByDiagonalMatrix(const QVector<double> &d)
		{
			Q_ASSERT(d.size() == rows);
			YaleSparseMatrix r(A, IA, JA);
			for (int i = 0; i < rows; ++i)
			{
				int endRow = IA[i + 1];
				for (int row = IA[i]; row < endRow; ++row)
				{
					int j = JA[row];
					r.A[row] *= d[i] * d[j];
				}
			}
			return r;
		}

		//not used, so not checked!
		double* GetDiagonal()
		{
			double* d = new double[rows];
			for (int i = 0; i < rows; ++i)
			{
				int endRow = IA[i + 1];
				for (int row = IA[i]; row < endRow; ++row)
				{
					int j = JA[row];
					if (j == i)
					{
						d[i] = A[row];
					}
				}
			}
			//Q_ASSERT(d.Length == RowCount);
			return d;
		}

		//not used, so not checked!
		YaleSparseMatrix PointWiseMultiply(YaleSparseMatrix B)
		{
			YaleSparseMatrix result(A,IA,JA);
			for (int i = 0; i < A.size(); ++i)
			{
				result.A[i] *= B.A[i];
			}
			return result;
		}

		bool Equals(YaleSparseMatrix other)
		{
			return A==other.A && IA==other.IA	&& JA == other.JA;
		}
};
}

#endif // YALESPARSEMATRIX_H
