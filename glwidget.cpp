#include "glwidget.h"
#include <GL/freeglut.h>

#define TOUCH true
#define MOUSE true
#define MOUSETOUCHID -1

#define MOUSEDIST 12
#define TOUCHDIST 35

enum {
	KEY_LEFT,
	KEY_RIGHT,
	KEY_FWD,
	KEY_BACK,
	KEY_UP,
	KEY_DOWN,
	KEY_SHIFT
};

GLWidget::GLWidget(QWidget *parent) :
	QGLWidget(parent),
	hasProblemSet(false),
	problemSet(nullptr),
	isolineRenderer(nullptr),
	frame_timer(new QTimer(this)),
	offset(),
	scale(1),
	width(-1),
	height(-1)
{
	setAttribute(Qt::WA_AcceptTouchEvents);
	makeCurrent(); // Ensure that openGL context is initialized.
	int a=1;
	glutInit(&a,nullptr);
	// Set a timer for updating the view periodically
	connect(frame_timer, SIGNAL(timeout()), SLOT(update_frame()));
	frame_timer->setSingleShot(true);
	scheduleRepaint();
	solveTimer.start();



	// Make the GL area focusable
	setFocusPolicy(Qt::StrongFocus);

	// Init initial key states
	for (size_t i = 0; i < sizeof(pressed) / sizeof(pressed[0]); i++)
		pressed[i] = false;

	QPixmap font("../normal.png");

	text_img = bindTexture(font, GL_TEXTURE_2D);
}

GLWidget::~GLWidget() {
}

void GLWidget::initializeGL() {
	QGLWidget::initializeGL();
#ifdef SPAMDEBUGMESSAGES
	qDebug("Version: %s", glGetString(GL_VERSION));
	qDebug("Vendor: %s", glGetString(GL_VENDOR));
	qDebug("Renderer: %s", glGetString(GL_RENDERER));
	qDebug("Extensions: %s", glGetString(GL_EXTENSIONS));
	GLint maxbuffers;
	glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS, &maxbuffers);
	qDebug("Max Color Attachments: %d", maxbuffers);
	Q_ASSERT(maxbuffers>=2);
#endif

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_SMOOTH);
	//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE); // GL_DECAL
	glDisable(GL_DEPTH_TEST);


	glClearColor(1, 1, 1, 0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GLWidget::viewport()
{
	glLoadIdentity();
	glScalef(scale,scale,1);
	glTranslatef(offset.rx(),-offset.ry(),0);
	scheduleRepaint();
}

void GLWidget::resizeGL(int w, int h) {
	QGLWidget::resizeGL(w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h); // update context viewport size
	glOrtho(0, w, 0, h, -1, 1); // set origin to left bottom corner.
	glMatrixMode(GL_MODELVIEW);
	viewport();
	width = w;
	height = h;
}

void GLWidget::reset()
{
	offset=QPointF();
	scale=1;
	viewport();
}



/**
 * Paints the contents of the OpenGL context.
 */
void GLWidget::paintGL() {
	QGLWidget::paintGL();
	glClear(GL_COLOR_BUFFER_BIT);

	if(problemSet==nullptr){return;}
	glColor3f(1,1,1);
	problemSet->LockMutex(QString::number(__LINE__));
	int red=0;
	int green=0;
	int blue=0;
	if(problemSet->hasTexture && true){
		glBindTexture(GL_TEXTURE_2D,problemSet->glTextureId);
		glEnable(GL_TEXTURE_2D);

		glBegin(GL_TRIANGLES);
#define V(x) problemSet->X.VertexData(x)
#define T(x) problemSet->triangles[x]
		for(TriangleIterator it=problemSet->mesh->Begin();it!=problemSet->mesh->End();it.Next()){
			if(!it.GetTriangleP().IsClockwise()){
				it.Triangle().color=TRED;
			}else{
				if(it.Triangle().color==TRED){it.Triangle().color=TNOCOLOR;}
			}
			if(it.Triangle().color==TGREEN){glColor3f(0,1,0);green++;}
			if(it.Triangle().color==TRED){glColor3f(1,0,0);red++;}
			if(it.Triangle().color==TBLUE){glColor3f(0,0,1);blue++;}
			if(!params.colors || it.Triangle().color==TNOCOLOR){glColor3f(1,1,1);}
			for(int j=0;j<3;j++){
				//double c=isolineRenderer->GetScalarValue(V(T(i)(j)).toQVector2D())/11.4*0.7+0.3;
				//glColor3f(c,c,c);
				glTexCoord2f(it[j].tx,it[j].ty);
				glVertex2f(it[j].x(),it[j].y());
			}
			Q_ASSERT(glGetError()==GL_NO_ERROR);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);

	if(params.smoothPoints){
		glEnable(GL_POINT_SMOOTH);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
	}
	if(params.showEdges){
		// Draw lines
		glColor4f(0,0,0,0.3);
		glLineWidth(params.lineWidth);
		glBegin(GL_LINES);
		TriangleIterator it=problemSet->mesh->Begin();
		for(;it!=problemSet->mesh->End();it.Next()){
			for(int i=0;i<3;i++){
				if(it.Triangle().stressEdge!=i || !params.colors){
					glVertex2d(it.GetEdge(i).first.x(),it.GetEdge(i).first.y());
					glVertex2d(it.GetEdge(i).second.x(),it.GetEdge(i).second.y());
				}
			}
		}
		glEnd();
	}
	if(params.showEdges && params.colors){
		// Draw lines
		glColor4f(1,0,0,1);
		glLineWidth(params.lineWidth+3);
		glBegin(GL_LINES);
		TriangleIterator it=problemSet->mesh->Begin();
		for(;it!=problemSet->mesh->End();it.Next()){
			for(int i=0;i<3;i++){
				if(it.Triangle().stressEdge==i){
					glVertex2d(it.GetEdge(i).first.x(),it.GetEdge(i).first.y());
					glVertex2d(it.GetEdge(i).second.x(),it.GetEdge(i).second.y());
				}
			}
		}
		glEnd();
	}
	if (problemSet->size() > 0 && params.showPoints) {
		// Draw points
		glColor3f(1,0,0);
		glPointSize(params.pointSize);
		glBegin(GL_POINTS);
		Vertex* v=problemSet->mesh->GetVertexPointer();
		for(;v!=nullptr;v=(Vertex*)v->next){
			if(!v->IsReal()){
				glVertex2d(v->x(),v->y());
			}
		}
		glEnd();
	}
	if(params.bridgesRidges){
		BridgesAndRidges(params.bridgeFactor);
	}

	if(params.showRealPoints){
		// Draw "real" points
		glColor3f(1,0,1);
		glPointSize(params.pointSize+5);
		glBegin(GL_POINTS);
		Vertex* v=problemSet->mesh->GetVertexPointer();
		for(;v!=nullptr;v=(Vertex*)v->next){
			if(v->IsReal()){
				glVertex2d(v->x(),v->y());
			}
		}
		glEnd();
		v=problemSet->mesh->GetVertexPointer();
		for(;v!=nullptr;v=(Vertex*)v->next){
			if(v->IsReal()){
				glColor4f(0,0,0,1);
				glRasterPos2f(v->x()+((params.pointSize+5)/2+2)*(1/scale),v->y());
				glutBitmapString(GLUT_BITMAP_HELVETICA_18, v->GetName());
			}
		}
	}
	if(params.showSelectedPoints){
		glBegin(GL_POINTS);
		for(int i=0;i<problemSet->TouchLocks.size();i++){ //draw selected/locked point(s)
			// Draw points
			if(problemSet->TouchLocks[i].first < -1){//normal lock
				glColor3f(0,0,1);
				glPointSize(params.pointSize+3);
			}
			else{
				glColor3f(0,1,0);
				glPointSize(params.pointSize);
			}
			Vertex* p=problemSet->TouchLocks[i].second.first;
			glVertex2d(p->x(),p->y());
		}
		glEnd();
	}
	//Draw error vecs
	if(params.arrowWeight>1e-6){
		glColor3f(1,0,0);
		glLineWidth(params.pointSize+3);
	}
	if(params.showIso){
		glColor3f(0,0,0);
		glLineWidth(2);
		//Draw isolines
		if(problemSet->hasAssociatedValues){
			isolineRenderer->DrawIsolines();
		}
	}
	DEBUG << "red: "<<red<< "   green: "<<green<<"    blue: "<<blue;
	problemSet->UnlockMutex(QString::number(__LINE__));
}

void GLWidget::BridgesAndRidges(double factor){
	if(!realVertsCalculated){
		realVerts.clear();
		auto it=problemSet->mesh->VBegin();
		for(;it!=problemSet->mesh->VEnd();it.Next()){
			if(it.Vertex().IsReal()){
				realVerts.append(it.VertexPointer());
			}
		}
		realVertsCalculated=true;
	}
	QVector<QVector2D> markers;
	for(int i=0;i<realVerts.size();i++){
		for(int j=i+1;j<realVerts.size();j++){
			bool faded=true;
			int ii=i,jj=j;
			if(problemSet->HasLockByPoint(realVerts[j])!=NOLOCK){faded=false;}
			if(problemSet->HasLockByPoint(realVerts[i])!=NOLOCK){faded=false;jj=i;ii=j;}
			double d=((*realVerts[i])-(*realVerts[j])).length();
			double dis=problemSet->mesh->Dis(realVerts[i]->ID(),realVerts[j]->ID());//*algoParams->realWeight;
			if(d<dis*(1-factor)){DrawRidge(*realVerts[ii],*realVerts[jj],1-d/dis,dis,faded,markers);}
			if(dis*(1+factor)<d){DrawBridge(*realVerts[ii],*realVerts[jj],1-dis/d,dis,faded,markers);}
		}
	}
	glColor3f(0,0,0);
	glPointSize(8);
	glBegin(GL_POINTS);
	for(int i=0;i<markers.size();i++){
		glVertex2d(markers[i].x(),markers[i].y());
	}
	glEnd();
}

void GLWidget::DrawRealDistLine(QVector2D from, QVector2D dir, double distance,QVector<QVector2D> &markers){
	glColor4f(1,1,1,0.8);
	glLineWidth(6);
	glBegin(GL_LINES);
	glVertex2d(from.x(),from.y());
	glVertex2d(dir.x(),dir.y());
	glEnd();
	glColor4f(0,0,0,1);
	if((from-dir).length()<distance){dir=(dir-from).normalized()*distance+from;}
	glLineWidth(1);
	glBegin(GL_LINES);
	DashedLine(from,dir);
	glEnd();
	QVector2D p=(dir-from).normalized()*distance+from;
	markers.append(p);
}

void GLWidget::DrawBridge(QVector2D a, QVector2D b,double d,double dis,bool faded, QVector<QVector2D> &markers){
	QVector2D mid=(a+b)/2;
	QVector2D offset=(b-a)*(d/4);
	QVector2D aa=mid+offset;
	QVector2D bb=mid-offset;
	float alpha=(faded)?0.4:1;
	if(!faded){DrawRealDistLine(a,b,dis,markers);}
	glColor4f(0.2,0.7,0.2,alpha);
	glLineWidth((faded)?3:4);
	glBegin(GL_LINES);
	glVertex2d(aa.x(),aa.y());
	glVertex2d(bb.x(),bb.y());
	glEnd();
}
void GLWidget::DrawRidge(QVector2D a, QVector2D b,double d,double dis,bool faded, QVector<QVector2D> &markers){
	QVector2D mid=(a+b)/2;
	QVector2D offset=(b-a)*(d/8);
	double tmp=offset.x();
	offset.setX(offset.y()*-1);
	offset.setY(tmp);
	glColor3f(1,0,0);
	QVector2D aa=mid+offset;
	QVector2D bb=mid-offset;
	float alpha=(faded)?0.4:1;
	if(!faded){
		DrawRealDistLine(a,b,dis,markers);
	}
	glLineWidth((faded)?3:4);
	glColor4f(0.8,0.2,0.2,alpha);
	glBegin(GL_LINES);
	glVertex2d(aa.x(),aa.y());
	glVertex2d(bb.x(),bb.y());
	glEnd();

}

void GLWidget::scheduleRepaint()
{
	if (!frame_timer->isActive())
		frame_timer->start(30);
}

void GLWidget::DrawVector(QVector2D vec,QVector2D origin){
	glVertex2d(origin.x(),origin.y());
	glVertex2d(origin.x()+vec.x(), origin.y()+vec.y());
}

void GLWidget::UpdateTouchList(const QList<QTouchEvent::TouchPoint>& points){
	problemSet->LockMutex(QString::number(__LINE__));
	//remove locks for lifted fingers
	for(int i=0;i<problemSet->TouchLocks.size();i++){
		bool found=false;
		for(int j=0;j<points.size();j++){
			if(points[j].id()==problemSet->TouchLocks[i].first){found=true;break;}
		}
		if(!found){//this touchId is no longer active, remove lock;
			problemSet->TouchLocks.removeAt(i);
		}
	}
	for (int i = 0;i<points.size(); i++) {
		//now we add or update locks
		Vertex* lock=problemSet->HasLock(points[i].id());
		if(lock != nullptr){//update position of selected point
			QPointF d=(points[i].pos() - points[i].lastPos())/scale;
			QVector2D v=QVector2D(lock->x()+d.x(),lock->y()-d.y());

			if(algoParams->guardedMovement){
				lock->MoveTo(v,algoParams->minProjectionHeight*algoParams->softHardFactor);
			}
			else{
				lock->setX(v.x());
				lock->setY(v.y());
			}
			//update lock
			problemSet->SetLock(points[i].id(),lock);
		}
		else{
			//new lock, find close point if it exists
			QPointF p = points[i].pos();
			p.setY(p.y()-height);
			p = p/scale-offset;

			double r=TOUCHDIST/scale;
			QVector2D v(p.x(),-p.y());
			VertexIterator it=problemSet->mesh->VBegin();
			for (; it!=problemSet->mesh->VEnd(); it.Next()) {
				double d=(v-it.Vertex()).length();
				if (d<r) {
					d=r;
					//if already locked, ignore
					int lock = problemSet->HasLockByPoint(it.VertexPointer());
					if(lock>=0){break;} //we can't move other touch locks untill we lift the other finger
					problemSet->SetLock(points[i].id(),it.VertexPointer());
					qDebug("Touching point (%f,%f) with touch Id %i",it.x(),it.y(),points[i].id());
					DEBUG << problemSet->TouchLocks;
					break;
				}
			}
		}
	}
	problemSet->UnlockMutex(QString::number(__LINE__));
	if(algoParams->solvesPerSec!= 0 && solveTimer.elapsed()<(1000/algoParams->solvesPerSec)) return;
	solveTimer.restart();
	if(params.autoSolve){
		emit(doSolve());}
	viewport();
}

bool GLWidget::event(QEvent * e){
	if(problemSet==nullptr){return QWidget::event(e);}
	switch(e->type()){
		case QEvent::TouchBegin:
		{
			DEBUG << "TouchBegin";
			const QList<QTouchEvent::TouchPoint>& points = static_cast<QTouchEvent*>(e)->touchPoints();
			UpdateTouchList(points);
			break;
		}
		case QEvent::TouchUpdate:
		{
			//DEBUG << "TouchUpdate";
			const QList<QTouchEvent::TouchPoint>& points = static_cast<QTouchEvent*>(e)->touchPoints();
			UpdateTouchList(points);
			//DEBUG << "Points: " << problemSet->TouchLocks;
			break;
		}
		case QEvent::TouchEnd:
		{
			DEBUG << "TouchEnd";
			//const QList<QTouchEvent::TouchPoint>& points = static_cast<QTouchEvent*>(e)->touchPoints();
			QList<QTouchEvent::TouchPoint> points; //we pass an empty list because for some reason the touch list still contains the last point after the touch end event
			UpdateTouchList(points);
			break;
		}
		default:
			return QWidget::event(e);
	}
	return true;
}

void GLWidget::mousePressEvent(QMouseEvent *event) {
	if(!MOUSE){return;}
	QWidget::mousePressEvent(event);
	if (event->button() == Qt::RightButton) {
		drag = event->pos();
		setCursor(QCursor(Qt::ClosedHandCursor));
		problemSet->Unlock(MOUSETOUCHID);
	}
	if (event->buttons() == Qt::MiddleButton) { //toggle lock
		QPointF p = event->pos();
		dragRight=p;
		p.setY(p.y()-height);
		p = p/scale-offset;

		double r=10/scale;
		QVector2D v(p.x(),-p.y());
		VertexIterator it=problemSet->mesh->VBegin();
		for (; it!=problemSet->mesh->VEnd(); it.Next()) {
			double d=(v-it.Vertex()).length();
			if (d<r) {
				d=r;
				if(problemSet->UnlockPoint(it.VertexPointer())){
					qDebug("UnLocked point %i (%f,%f)\n",it.Vertex().ID(),it.x(),it.y());
				}else{
					problemSet->SetLock(-it.Vertex().ID()-2,it.VertexPointer());
					qDebug("Locked point %i (%f,%f)\n",it.Vertex().ID(),it.x(),it.y());
				}
				break;
			}
		}
	}
	if (event->buttons() == Qt::LeftButton) {
		if(problemSet==nullptr){return;}
		QPointF p = event->pos();
		dragRight=p;
		p.setY(p.y()-height);
		p = p/scale-offset;
		DEBUG << "Click pos = "<<QPointF(p.x(),-p.y())<<isolineRenderer->GetScalarValue(QVector2D(p.x(),-p.y()));
		double r=MOUSEDIST/scale;
		QVector2D v(p.x(),-p.y());
		VertexIterator it=problemSet->mesh->VBegin();
		for (; it!=problemSet->mesh->VEnd(); it.Next()) {
			double d=(v-it.Vertex()).length();
			if (d<r) {
				d=r;
				//remove old lock if needed
				problemSet->Unlock(MOUSETOUCHID);
				problemSet->SetLock(MOUSETOUCHID,it.VertexPointer());

				qDebug("Selected point %i (%f,%f) with mouse\n",it.Vertex().ID(),it.x(),it.y());
				DEBUG << "VALUE: "<<it.Vertex().value;
				break;
			}
		}
		setCursor(QCursor(Qt::ClosedHandCursor));
	}
	viewport();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event) {
	if(!MOUSE){return;}
	QWidget::mouseReleaseEvent(event);
	if (event->buttons() == 0) {
		setCursor(QCursor(Qt::ArrowCursor));
	}
	viewport();
}

/** Uses mouse input to do (...).
 * The following actions are assigned to mouse buttons:
 *  - Left: drag viewport
 *  - Middle: -
 *  - Right: select vertex
 * The following modifiers are used:
 *  - Shift: -
 *  - Control: -
 */
void GLWidget::mouseMoveEvent(QMouseEvent *event) {
	if(problemSet==nullptr){return;}
	if(!MOUSE){return;}

	if (event->modifiers() & Qt::ShiftModifier) {}
	if (event->modifiers() & Qt::ControlModifier) {}

	if (event->buttons() == Qt::RightButton) {
		QWidget::mouseMoveEvent(event);
		QPointF d = (event->pos() - drag); // Compute moved distance
		offset+=d/scale;
		viewport();
		drag = event->pos();
	}
	else if (event->buttons() == Qt::MiddleButton) {}
	else if (event->buttons() == Qt::LeftButton) {
		if(algoParams->solvesPerSec!= 0 && solveTimer.elapsed()<(1000/algoParams->solvesPerSec)) return;
		solveTimer.restart();
		QWidget::mouseMoveEvent(event);
		problemSet->LockMutex(QString::number(__LINE__));
		QPointF d = (event->pos() - dragRight)/scale;
		Vertex* lock=problemSet->HasLock(MOUSETOUCHID);
		if(lock!=nullptr){ //we are moving a point
			QVector2D v=QVector2D(lock->x()+d.x(),lock->y()-d.y());

			if(algoParams->guardedMovement){
				lock->MoveTo(v,algoParams->minProjectionHeight*algoParams->softHardFactor);
			}
			else{
				lock->setX(v.x());
				lock->setY(v.y());
			}
			//update lock
			problemSet->SetLock(MOUSETOUCHID,lock);
		}
		problemSet->UnlockMutex(QString::number(__LINE__));
		if(params.autoSolve){
			emit(doSolve());}
		viewport();

		dragRight=event->pos();
	}
}

/** Uses mouse wheel scrolling for zooming in and out.
 */
void GLWidget::wheelEvent(QWheelEvent *event) {
	QWidget::wheelEvent(event);
	QPointF p = event->pos();
	p.ry()-=height;
	offset-=p/scale;
	scale *= pow(2,event->delta()/720.);
	offset+=p/scale;
	viewport();
}

/** Handles the key event to keep track of button press state.
 * Returns true if pressing the key requires periodic updates of the
 * OpenGL context.
 */
bool GLWidget::setKey(int keycode, bool state) {
	switch(keycode) {
		case Qt::Key_R:
			reset();
			return true;
		case Qt::Key_Left:
			pressed[KEY_LEFT] = state;
			return true;
		case Qt::Key_Right:
			pressed[KEY_RIGHT] = state;
			return true;
		case Qt::Key_A:
			pressed[KEY_FWD] = state;
			return true;
		case Qt::Key_Z:
			pressed[KEY_BACK] = state;
			return true;
		case Qt::Key_Down:
			pressed[KEY_DOWN] = state;
			return true;
		case Qt::Key_Up:
			pressed[KEY_UP] = state;
			return true;
		case Qt::Key_Shift:
			pressed[KEY_SHIFT] = state;
			return false;
		default:
			return false;
	}
}

void GLWidget::keyPressEvent(QKeyEvent *e) {
	if (!e->isAutoRepeat() && setKey(e->key(), true)) {
		scheduleRepaint();
	} else {
		QWidget::keyPressEvent(e);
	}
}

void GLWidget::keyReleaseEvent(QKeyEvent *e) {
	if (!e->isAutoRepeat() && setKey(e->key(), false)) {
	} else {
		QWidget::keyReleaseEvent(e);
	}
}

void GLWidget::update_frame() {
	bool v=false;
	float s=25.f/scale;
	offset.rx()-=width/2.f/scale;
	offset.ry()+=height/2.f/scale;
	if (pressed[KEY_FWD]) {
		scale*=1.1;
		v=true;
	}
	if (pressed[KEY_BACK]) {
		scale/=1.1;
		v=true;
	}
	offset.rx()+=width/2.f/scale;
	offset.ry()-=height/2.f/scale;
	if (pressed[KEY_LEFT]) {
		offset.rx()+=s;
		v=true;
	}
	if (pressed[KEY_RIGHT]) {
		offset.rx()-=s;
		v=true;
	}
	if (pressed[KEY_DOWN]) {
		offset.ry()-=s;
		v=true;
	}
	if (pressed[KEY_UP]) {
		offset.ry()+=s;
		v=true;
	}
	if (v) viewport();
	updateGL(); // only repaint call.
}

void GLWidget::setVisParams(VisParams & _params){
	params=_params	;
	if(isolineRenderer!=nullptr){
		isolineRenderer->gridSize=params.isoGrid;}
	scheduleRepaint();
}

void GLWidget::setProblemSet(ProblemSet *p, AlgoParams *a) {
	problemSet = p;
	algoParams=a;
	hasProblemSet=true;
	isolineRenderer=new DisIsolines(p);
	realVertsCalculated=false;
	scheduleRepaint();
}

