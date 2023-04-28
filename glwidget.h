#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QApplication>
#include "structures.h"
#include "isolines.h"

#include <cmath>
#include <cassert>

#include <QTime>
#include <QTimer>
#include <QDebug>
#include <QSettings>
#include <QMouseEvent>
#include <limits>


#if defined(_WIN32) || defined(_WIN64)
# define GL_MAX_COLOR_ATTACHMENTS          0x8CDF
#else
#include <QtOpenGL>
#endif

class QTimer;

class GLWidget : public QGLWidget {
		Q_OBJECT
	public:
		explicit GLWidget(QWidget *parent = 0);
		virtual ~GLWidget();
		VisParams params;
		AlgoParams * algoParams;
		bool hasProblemSet;
		ProblemSet * problemSet;
		void setProblemSet(ProblemSet * p, AlgoParams *a);
	protected:
		virtual void initializeGL();
		virtual void resizeGL(int w, int h);
		virtual void paintGL();
		virtual void mousePressEvent(QMouseEvent *e);
		virtual void mouseReleaseEvent(QMouseEvent *e);
		virtual void mouseMoveEvent(QMouseEvent *e);
		virtual void keyPressEvent(QKeyEvent *e);
		virtual void keyReleaseEvent(QKeyEvent *e);
		virtual void wheelEvent(QWheelEvent *e);
		virtual bool event(QEvent *e); // for touch events
	private:
		DisIsolines * isolineRenderer;
		QTime solveTimer;
		QTimer * frame_timer;
		bool pressed[7]; // keys that are pressed.
		bool setKey(int keycode, bool state);
		QPointF drag, dragRight;
		QPointF offset;
		float scale;

		int width;
		int height;
		int text_img;
		void viewport();

		void UpdateTouchList(const QList<QTouchEvent::TouchPoint>& points);
		void DrawVector(QVector2D vec, QVector2D origin=QVector2D(0,0));
		void BridgesAndRidges(double factor);
		void DrawRidge(QVector2D a, QVector2D b, double d, double dis, bool faded, QVector<QVector2D> &markers);
		void DrawBridge(QVector2D a, QVector2D b, double d, double dis, bool faded, QVector<QVector2D> &markers);
		void DashedLine(QVector2D a, QVector2D b, double dashlen=2){
			int dashes=int((a-b).length()/dashlen)/2;
			QVector2D dashvec=(b-a).normalized()*dashlen;
			double offset=(a-b).length()-(dashlen*dashes*2-1);
			a+=(b-a).normalized()*offset;
			for(int i=0;i<dashes;i++){
				glVertex2d(a.x(),a.y());
				a+=dashvec;
				glVertex2d(a.x(),a.y());
				a+=dashvec;
			}
		}
		void DrawRealDistLine(QVector2D from, QVector2D direction, double distance, QVector<QVector2D> &markers);
		QVector<Vertex*> realVerts;
		bool realVertsCalculated=false;
	private slots:
		void update_frame();
	public slots:
		void scheduleRepaint();
		void reset();
		void setVisParams(VisParams & params);
	signals:
		void changeRunning(bool);
		void doSolve();
};


#endif // GLWIDGET_H
