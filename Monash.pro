#-------------------------------------------------
#
# Project created by QtCreator 2013-11-11T01:45:17
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Monash
TEMPLATE = app
LIBS +=-lglut

SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp \
    config.cpp \
    solver.cpp \
    generator.cpp \
    worker.cpp \
    timing.cpp \
    random.cpp \
    SparseMajorization.cpp \
    del_impl.cpp \
    structures.cpp \
    londonbikes.cpp \
    isolines.cpp \
    dynamicmesh.cpp \
    geom.cpp \
    DynamicMesh_impl.cpp

HEADERS  += mainwindow.h \
    structures.h \
    glwidget.h \
    config.h \
    solver.h \
    generator.h \
    worker.h \
    timing.h \
    random.h \
    YaleSparseMatrix.h \
    SparseMajorization.h \
    HundredNodeTriangulation.h \
    ProjectionMethods.h \
    Texture.h \
    del_interface.h \
    triangle.h \
    triangle_impl.h \
    del_interface.hpp \
    dpoint.h \
    isolines.h \
    dynamicmesh.h \
    DynamicMesh_impl.h \
    geom.h

QMAKE_CXXFLAGS += -std=c++11

OTHER_FILES +=
