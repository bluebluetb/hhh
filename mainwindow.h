#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#include <QMainWindow>
#include <QLabel>
#include <QSettings>
#include <QDebug>
#include <QTimer>
#include <QThread>
#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QDockWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QPushButton>

#include <sstream>

#include "glwidget.h"
#include "worker.h"
#include "generator.h"
#include "config.h"
#include "solver.h"
#include "timing.h"
#include "config.h"
#include "structures.h"
#include "SparseMajorization.h"

class Task;
class QThread;
class OptionItem;
class Config;
class GLWidget;

class Window : public QMainWindow {
				Q_OBJECT
public:
	Window();
		~Window();
	virtual void closeEvent(QCloseEvent *e);

private:
	Config *createConfig(QString name);
	void createMenu();
	void createConfigPanels();
	void createMenuOptions();
	void createContext();
		void loadSettings();
	FileChooserItem * imageSelect=nullptr;
	QLabel * info;
//	DoubleItem *visParamItem, *visParamItem2;
	struct {
		QMenu *file;
		QMenu *window;
		QMenu *help;
	} menu;
	ButtonItem * go;
	GLWidget *gl;
	QThread *workThread;
	QTimer continuousSolveTimer;
private slots:
	void aboutQt();
	void generateProblem();
	void solve();
	void createWorker();
	void VisParamChanged();
	void ProcessTaskInfo(TaskInfo*);
	void SubDiv();
	void Fixup();
	void Fixup2();
	void DeleteSelected();
signals:
		void solved(TaskInfo*);
		void newTask(Task*);
		void quitWorkerThread();
		void stopWorkerTask();
};

#endif // FRAMEWORK_H
