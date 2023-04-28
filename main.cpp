#include <QApplication>
#include "mainwindow.h"

int TaskInfo::solves=0;

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	Window w;
	w.show();
	return a.exec();
}
