#ifndef WORKER_H
#define WORKER_H
#include <QObject>
#include <QHash>
#include "structures.h"

class Task : public QObject, public BaseTask {
		Q_OBJECT
	public:
		Task(ProblemSet &p,  AlgoParams &algoParams);
		~Task(){}
		void create_edge(Edge e);
	public slots:
		void do_quit();
	signals:
		void working(bool);
		void add_edge(Edge e);
		friend class Worker;
};

class Worker : public QObject {
		Q_OBJECT
	public:
		explicit Worker();
		virtual ~Worker();
	public slots:
		void solve(Task *T);
	signals:
		void solved(TaskInfo*);
};

#endif // WORKER_H
