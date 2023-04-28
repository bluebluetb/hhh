#ifndef ISOLINES_H
#define ISOLINES_H

#include<QPair>
#include<QRectF>
#include<GL/gl.h>
#include"structures.h"
#include"del_interface.h"

class Isolines
{
	private:
		virtual double GetScalarValue(QVector2D pos)=0;
		void DrawIsoline(float value);
	public:
		QRectF window;
		int gridSize;
		int nLines;
		Isolines(QRectF window, int gridSize, int nLines=5);
		virtual void DrawIsolines();

};

class DisIsolines : public Isolines{
	private:
		ProblemSet * p;
		QVector2D topLeft, bottomRight;
		void UpdateVals();
		QVector<DMesh::Triangle> triangulation;
		QVector<Vertex*> realPoints;
	public:
		DisIsolines(ProblemSet *p):Isolines(QRectF(),30),p(p){}
		virtual double GetScalarValue(QVector2D pos);
		virtual void DrawIsolines();
};

#endif // ISOLINES_H
