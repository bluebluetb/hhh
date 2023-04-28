#include <QThread>
#include "worker.h"
#include "solver.h"
#include "timing.h"


Task::Task(ProblemSet &p, AlgoParams &algoParams)
		: QObject(), BaseTask(p,algoParams){}

/**
 * Creates a new edge.
 */
void Task::create_edge(Edge e) {
    emit add_edge(e);
		return;
}


Worker::Worker(): QObject() {
    qDebug("Worker [%p] Created", this);
}

Worker::~Worker() {
    qDebug("Worker [%p] Deleted", this);
}

void Worker::solve(Task *T) {
		Timing::reset();
		//qDebug("Starting solver for %d points using %s",T->p.size(),Solver::solvers[T->algoParams.algoID].name);
    emit T->working(true);
		Solver::solvers[T->algoParams.algoID].function(T);
		//qDebug("Solver finished");
		T->info.executionTime = Timing::elapsed();
		//qDebug("Time elapsed: %.15lfms",T->info.executionTime);
		T->info.solves++;
		emit T->working(false);
		emit solved(new TaskInfo(T->info));//create new TaskInfo to make sure it isnt destroyed before handling info

//    T->deleteLater();
		delete T;
}

void Task::do_quit() {
		p.quit=true;
}
