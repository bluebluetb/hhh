
#include "mainwindow.h"




namespace {

ProblemSet p;

OutputItem * debugOut;
VisParams params=VisParams(true);
GenParams genParams=GenParams(0,"../data/miniscale_std_with_grid_r16_20pc.png");
//GenParams genParams=GenParams(0,"/home/q/Monash/australia.jpg");
AlgoParams algoParams=AlgoParams(1);
}

Window::Window() {
	setDockOptions(AllowNestedDocks | AllowTabbedDocks);
	qRegisterMetaType<ProblemSet>("ProblemSet");
	qRegisterMetaType<QTextCursor>("QTextCursor");
	qRegisterMetaType<QTextCharFormat>("QTextCharFormat");
	qRegisterMetaType<Vertex>("Vertex");
	qRegisterMetaType<Edge>("Edge");
	qRegisterMetaType<Task*>("Task*");



	createContext();
	createWorker();
	createMenu();
	createConfigPanels();
	createMenuOptions();
	loadSettings();
	VisParamChanged();
}

Window::~Window() {
	// Let the worker exit its eventloop and wait for it to terminate
	emit quitWorkerThread();
	workThread->exit();
}

void Window::closeEvent(QCloseEvent *e) {
	QSettings s;
	s.setValue("geometry", saveGeometry());
	s.setValue("windowState", saveState());
	QWidget::closeEvent(e);
}

void Window::loadSettings() {
	resize(640, 480);
	QSettings s;
	restoreGeometry(s.value("geometry").toByteArray());
	restoreState(s.value("windowState").toByteArray());
}



void myMessageOutput(QtMsgType type, const char *msg)
{
	//in this function, you can write the message to any stream!
	char buffer [1024];
	sprintf(buffer, "%8d",int(Timing::elapsed()+0.5));
	QString tmp=buffer;
	tmp.append(" > ");
	tmp.append(msg);

	switch (type) {
		case QtDebugMsg:
			fprintf(stderr, "Debug: %s\n", msg);
			break;
		case QtWarningMsg:
			fprintf(stderr, "Warning: %s\n", msg);
			break;
		case QtCriticalMsg:
			fprintf(stderr, "Critical: %s\n", msg);
			break;
		case QtFatalMsg:
			fprintf(stderr, "Fatal: %s\n", msg);
			abort();
	}
}

void Window::createWorker() {
	// Create a thread for our worker.
	workThread = new QThread(this);
	workThread->start();
	// Create the worker.
	Worker * worker = new Worker();
	// Push it to his thread (Note that the thread object itself lives on our thread)
	worker->moveToThread(workThread);
	// Create the message queues
	connect(this, SIGNAL(newTask(Task*)), worker, SLOT(solve(Task*)));
	connect(worker, SIGNAL(solved(TaskInfo*)), this, SIGNAL(solved(TaskInfo*)));
	connect(worker, SIGNAL(solved(TaskInfo*)), this, SLOT(ProcessTaskInfo(TaskInfo*)));

	// Setup normal thread quiting sequence.
	connect(this, SIGNAL(quitWorkerThread()), worker, SLOT(deleteLater()));
	connect(worker, SIGNAL(destroyed(QObject*)), workThread, SLOT(quit()));
	qDebug("Created new Worker Thread.");
}

void Window::generateProblem() {
	delete p.mesh;
	p.mesh=new DynamicMesh();
	p = Generator::generators[genParams.genID].function(genParams);
    gl->setProblemSet(&p,&algoParams);
	imageSelect->SetFile(genParams.image);
	TaskInfo::solves=0;
}

void Window::solve() {
	if(p.size()==0){generateProblem();}
	if(!algoParams.continuousSolve){
		go->setText("Solve");
		continuousSolveTimer.stop();
		go->enable();
	}
	else{
		if(!continuousSolveTimer.isActive()){
			continuousSolveTimer.start(1000/algoParams.solvesPerSec);
			go->setText("Solving...");
		}
		if(go->isEnabled()){
			go->disable();
		}
	}
	Task * T = new Task(p, algoParams);
	VisParamChanged();
	T->moveToThread(workThread);

	connect(this,SIGNAL(stopWorkerTask()),T,SLOT(do_quit()),Qt::DirectConnection);

	emit newTask(T);
}

void Window::SubDiv(){
	p.LockMutex();
	//Solver::SubdivideTriangles(p);
	Solver::SubdivAndMerge(algoParams,p);
	gl->scheduleRepaint();
	p.UnlockMutex();
}

void Window::Fixup(){
	p.LockMutex();
	Solver::FlipTriangles(p,algoParams,true);
	p.UnlockMutex();
}

void Window::Fixup2(){
	p.LockMutex();
	Solver::FlipTriangles(p,algoParams,false);
	p.UnlockMutex();
}

void Window::DeleteSelected(){
	p.LockMutex();
	for(int i=p.TouchLocks.size()-1;i>=0;i--){
		if(p.TouchLocks[i].first>=-1 && !p.TouchLocks[i].second.first->IsReal()){
			auto vIt=p.mesh->VIterator(p.TouchLocks[i].second.first);
			vIt.RemoveAndMerge();
			p.TouchLocks.removeAt(i);
		}
	}
	gl->scheduleRepaint();
	p.UnlockMutex();
}
void Window::createConfigPanels() {
	{
		Config * c = createConfig("Generator");
		OptionItem * list = new OptionItem(c, "Generator", &genParams.genID);
		for (int i=0; Generator::generators[i].name; i++) {
			list->addItem(Generator::generators[i].name);
		}
		IntegerItem * s = new IntegerItem(c, "Seed", &genParams.seed);
		s->setRange(0,10000);
		IntegerItem * m = new IntegerItem(c, "Vertices", &genParams.count);
		m->setRange(0,1000);
		m->setSingleStep(1);
		IntegerItem * g = new IntegerItem(c, "Grid size", &genParams.grid);
		g->setRange(0,1000);
		g->setSingleStep(1);
		/*IntegerItem * r = new IntegerItem(c, "Max coordinate", &genParams.maxc);
		r->setRange(0,10000);
		r->setSingleStep(50);*/
        imageSelect=new FileChooserItem(c, "Image", &genParams.image, "Image Files (*.png *.jpg *.gif *.bmp)", "../Monash");

		BoolItem * grid = new BoolItem(c, "Overlay grid", &genParams.overlayGrid);
		IntegerItem * gx = new IntegerItem(c, "Grid vert lines", &genParams.overlayGridX);
		gx->setRange(3,30);
		gx->setSingleStep(1);
		IntegerItem * gy = new IntegerItem(c, "Grid hor lines", &genParams.overlayGridY);
		gy->setRange(3,30);
		gy->setSingleStep(1);


		ButtonItem * go = new ButtonItem(c, "Generate vertices");
		connect(go,SIGNAL(clicked(bool)),this,SLOT(generateProblem()));
		c->done();
	}
	{
		Config * c = createConfig("Visualization");

		BoolItem * points = new BoolItem(c, "Show points", &params.showPoints);
		BoolItem * rpoints = new BoolItem(c, "Show real", &params.showRealPoints);

		BoolItem * lpoints = new BoolItem(c, "Show locked", &params.showLockedPoints);
		BoolItem * spoints = new BoolItem(c, "Show selected", &params.showSelectedPoints);

		//BoolItem * smooth = new BoolItem(c, "Smooth points", &params.smoothPoints);
		BoolItem * edges = new BoolItem(c, "Show Edges", &params.showEdges);
		BoolItem * iso = new BoolItem(c, "Show Isolines", &params.showIso);
		BoolItem * bridges=new BoolItem(c, "BridgesRidges", &params.bridgesRidges);
		DoubleItem * bridgef=new DoubleItem(c, "Bridge tollerance", &params.bridgeFactor);
		bridgef->setRange(0,1);
		bridgef->setDecimals(1);
		bridgef->setSingleStep(0.1);
		connect(bridgef,SIGNAL(valueChanged(double)),this,SLOT(VisParamChanged()));

		connect(points,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		connect(rpoints,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		connect(lpoints,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		connect(spoints,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		connect(bridges,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		//connect(smooth,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		connect(edges,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		connect(iso,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));


		IntegerItem * pointSize = new IntegerItem(c,"Point size", &params.pointSize);
		pointSize->setRange(1,20);
		pointSize->setValue(params.pointSize);
		connect(pointSize,SIGNAL(valueChanged(int)),this,SLOT(VisParamChanged()));

		IntegerItem * lineWidth = new IntegerItem(c,"Line width", &params.lineWidth);
		lineWidth->setRange(1,10);
		lineWidth->setValue(params.lineWidth);
		connect(lineWidth,SIGNAL(valueChanged(int)),this,SLOT(VisParamChanged()));

		IntegerItem * isoGrid = new IntegerItem(c,"isoGridSize", &params.isoGrid);
		isoGrid->setRange(10,500);
		isoGrid->setValue(params.isoGrid);
		connect(isoGrid,SIGNAL(valueChanged(int)),this,SLOT(VisParamChanged()));

		/*DoubleItem * a = new DoubleItem(c, "Arrow mult", &params.arrowWeight);
		a->setRange(0,20);
		a->setSingleStep(0.1);
		a->setDecimals(1);
		connect(a,SIGNAL(valueChanged(double)),this,SLOT(VisParamChanged()));*/

		BoolItem * colors = new BoolItem(c, "Colors", &params.colors);
		connect(colors,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));

		BoolItem * solve = new BoolItem(c, "Auto solve", &params.autoSolve);
		connect(solve,SIGNAL(valueChanged(bool)),this,SLOT(VisParamChanged()));
		c->done();
	}
	{
		Config * c = createConfig("Solver");
		OptionItem * list = new OptionItem(c, "Solver", &algoParams.algoID);
		for (int i=0; Solver::solvers[i].name; i++) {
			list->addItem(Solver::solvers[i].name);
		}
		IntegerItem * l = new IntegerItem(c, "Limit", &algoParams.limit);
		l->setRange(0,100);
		l->setSingleStep(1);
		IntegerItem * s = new IntegerItem(c, "Max solves/s", &algoParams.solvesPerSec);
		s->setRange(0,100);
		s->setSingleStep(1);
		DoubleItem * mh = new DoubleItem(c, "Min HeightF", &algoParams.minProjectionHeight);
		mh->setRange(0,1);
		mh->setDecimals(2);
		mh->setSingleStep(0.05);
		DoubleItem * maxh = new DoubleItem(c, "Max HeightF", &algoParams.maxProjectionHeight);
		maxh->setRange(1,5);
		maxh->setDecimals(1);
		maxh->setSingleStep(0.1);
		new BoolItem(c, "Single minHeight", &algoParams.singleMinHeight);
		DoubleItem * mm = new DoubleItem(c, "SoftHard Factor", &algoParams.softHardFactor);
		mm->setRange(0,1);
		mm->setDecimals(2);
		mm->setSingleStep(0.05);
		new BoolItem(c, "Complete project only", &algoParams.onlyCompleteProjections);
		new BoolItem(c, "Guarded movement", &algoParams.guardedMovement);
		BoolItem* cSolve=new BoolItem(c, "Continuos Solve", &algoParams.continuousSolve);
		connect(cSolve,SIGNAL(valueChanged(bool)),SLOT(solve()));

		DoubleItem * sw = new DoubleItem(c, "Extra selectionW", &p.selectionWeight);
		sw->setDecimals(1);
		sw->setRange(0,10);
		sw->setSingleStep(0.5);
		DoubleItem * w = new DoubleItem(c, "Real Weight", &algoParams.realWeight);
		w->setRange(0,99);
		w->setSingleStep(0.1);
		w->setDecimals(1);

		DoubleItem * ma = new DoubleItem(c, "Min Triangle Area", &p.minTriangleArea);
		ma->setRange(0,500);
		ma->setSingleStep(20);
		ma->setDecimals(0);



		go = new ButtonItem(c, "Solve");
		connect(go,SIGNAL(clicked(bool)),this,SLOT(solve()));
		connect(gl,SIGNAL(doSolve()),this,SLOT(solve()),Qt::QueuedConnection);//Qt::DirectConnection);
		continuousSolveTimer.setSingleShot(false);
		connect(&continuousSolveTimer,SIGNAL(timeout()),SLOT(solve()));
		connect(go,SIGNAL(clicked(bool)),go,SLOT(disable()));
		connect(this,SIGNAL(solved(TaskInfo*)),go,SLOT(enable()));
		connect(this,SIGNAL(stopWorkerTask()),go,SLOT(enable()));
		LabelItem * infoItem = new LabelItem(c);
		info=infoItem->getLabel();


		c->done();

		ButtonItem * kill = new ButtonItem(c, "Stop current computation");
		connect(kill,SIGNAL(clicked(bool)),this,SIGNAL(stopWorkerTask()));

	}
	{
		Config * c = createConfig("Fixup");
		ButtonItem * s = new ButtonItem(c, "SubDiv/Merge");
		connect(s,SIGNAL(clicked(bool)),this,SLOT(SubDiv()));
		ButtonItem * go = new ButtonItem(c, "Fixup");
		connect(go,SIGNAL(clicked(bool)),this,SLOT(Fixup()));
		ButtonItem * go2 = new ButtonItem(c, "Fixup unguarded");
		connect(go2,SIGNAL(clicked(bool)),this,SLOT(Fixup2()));
		ButtonItem * del = new ButtonItem(c,"Delete selected");
		connect(del,SIGNAL(clicked(bool)),this,SLOT(DeleteSelected()));
	}
	{
		Config * c = createConfig("Debug output");
		debugOut = new OutputItem(c);
		//qInstallMsgHandler(myMessageOutput);/////////////////////////////
	}
}

void Window::ProcessTaskInfo(TaskInfo *i){
	using namespace std;
	stringstream ss;
	ss<<"Time elapsed: "<<i->executionTime<<"\n";
	ss<<"Last stress: "<<i->lastStress<<"\n";
	ss<<"Total solves: "<<i->solves<<"\n";
	info->setText(QString::fromStdString(ss.str()));
	delete i; //i was created with new so it wouldn't get destroyed before the signal was processed
}

void Window::createMenu() {
	// Add menu options
	QMenuBar *bar = menuBar();
	menu.file   = bar->addMenu("File");
	menu.window = bar->addMenu("Window");
	menu.help   = bar->addMenu("Help");
}

void Window::createMenuOptions() {
	QAction *a;

	menu.file->addSeparator();
	a = new QAction("Quit", this);
	connect(a, SIGNAL(triggered()), SLOT(close()));
	menu.file->addAction(a);

	a = new QAction("About Qt4", this);
	connect(a, SIGNAL(triggered()), SLOT(aboutQt()));
	menu.help->addAction(a);
}

void Window::aboutQt() {
	QMessageBox::aboutQt(this);
}


Config *Window::createConfig(QString name) {
	// Create configuration toolbar
	QDockWidget *dock = new QDockWidget(name, this);
	dock->setObjectName(name);
	addDockWidget(Qt::LeftDockWidgetArea, dock);
	menu.window->addAction(dock->toggleViewAction());
	Config *c = new Config(dock);
	dock->setWidget(c);
	return c;
}

void Window::createContext() {
	// Create OpenGL context
	gl = new GLWidget(this);
	setCentralWidget(gl);
}

void Window::VisParamChanged(){
	gl->setVisParams(params);
	gl->scheduleRepaint();
}
