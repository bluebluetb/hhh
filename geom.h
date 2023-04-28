#ifndef GEOM_H
#define GEOM_H

#include <QVector2D>
#include <QVector3D>

#define EPS 1e-6

bool LinesParallel(QVector2D a, QVector2D b, QVector2D c, QVector2D d);
bool LinesCollinear(QVector2D a, QVector2D b, QVector2D c, QVector2D d);
bool SegmentsIntersect(QVector2D a, QVector2D b, QVector2D c, QVector2D d);
bool SegmentLineIntersect(QVector2D p1, QVector2D p2, QVector2D a, QVector2D b);
QVector2D ComputeLineIntersection(QVector2D a, QVector2D b, QVector2D c, QVector2D d);
bool SegmentLineIntersect(QVector2D p1, QVector2D p2, QVector2D a, QVector2D b);
QVector2D RotateCW90(QVector2D p);
QVector2D RotateCCW90(QVector2D p);

#endif // GEOM_H
