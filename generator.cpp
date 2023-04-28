#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <QDebug>

#include "random.h"
#include "generator.h"
#include "Texture.h"
#include "del_interface.h"
#include "quadtree_mesh.h"

#include <QFile>
#include <QTextStream>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

using namespace Generator;

#include "SparseMajorization.h"
#include "worker.h"
void AddSolutionAsMovement(ProblemSet &realp, ProblemSet &p){
	VertexContainer points(1);
	QVector<Vertex *> origs;
	for(auto it=realp.mesh->VBegin();it!=realp.mesh->VEnd();it.Next()){
		points.append(it.VertexPointer());
		origs.append(new Vertex(it.Vertex()));
	}
	AlgoParams params;
	Task * T = new Task(realp,params);
	T->algoParams.limit=0;
	T->algoParams.realWeight=1;
	T->algoParams.guardedMovement=false;
	TriangleProjection::SparseMajorization majorization(T,points);
	majorization.Run(T,points);
	//points now contains end positions for our original vertices
	for(int i=0;i<origs.size();i++){
		Vertex * v=origs[i];
		p.mesh->AddVertex(v);
		QVector2D target;target.setX(points(0,i));target.setY(points(1,i));
		p.moves.append(MovementData(v,(*v),target));
	}
	p.mesh->Triangulate();
}

ProblemSet DFile(GenParams & params){//this is going to be the ID of the first real data point.
	ProblemSet p;
	QVector<Vertex> points;
	srand(params.seed);
	int N=params.grid;
	int M=params.grid;
	QString fileName="../../data/in.txt";
	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
		qDebug() << "error opening file \"" << fileName << "\"!\n";;
		return p;
	}
	QTextStream in(&file);
	QString imageFile="../../data/"+in.readLine();
	QImage img;
	if(!img.load(imageFile)){
		qDebug() << "error loading image file" << imageFile <<"!\n";
		return p;
	}
	int w=img.width();
	int h=img.height();

	int V=0;
	//add data points
	while (!in.atEnd()) {
		QString line = in.readLine();
		if (line.startsWith("EDGES",Qt::CaseInsensitive)) break;
		QTextStream l(line.toStdString().c_str());
		double x,y;
		l >> x;
		l >> y;
		QString label;
		l >> label;
		Vertex v;
		v.SetReal(true);
		v.SetName(label);
		v.setX(x);
		v.setY(y);
		//v.value=c;
		points.append(v);
		V++;
	}

	//build grid
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D((double(i)/(M-1))*w,(double(j)/(N-1))*h));
		}
	}

	p.mesh->SetVertices(points);
	int offset=p.mesh->VBegin().Vertex().ID();
	//add "Real" dissimilarities
	while (!in.atEnd()) {
		int a, b;
		double v;
		QString line = in.readLine();
		QTextStream l(line.toStdString().c_str());
		l >> a >> b >> v;
		p.mesh->SetDis(a+offset,b+offset,v);

	}
	file.close();
	params.image=imageFile;
	TextureLoader::FromFile(params,p);
	return p;
}

/*
ProblemSet OneTriangle(GenParams & p){
	double ** D =new double*[3];
	for(int i=0;i<3;i++){D[i]=new double[3];}
	D[0][0]= 0;D[0][1]=10;D[0][2]=10;
	D[1][0]=10;D[1][1]= 0;D[1][2]=10;
	D[2][0]=10;D[2][1]=10;D[2][2]= 0;

	// random starting positions
	srand(time(0));
	Matrix2xN X;
	X.resize(3);
	for(unsigned int i=0;int(i)<X.rows();i++){
		X.Set(i,Vertex(rand()%p.maxc, rand()%p.maxc));
	}

	return SetEdges(ProblemSet(X,D));
}

ProblemSet LittleGraph(GenParams& p){
	double** D = new double*[5];
	for(int i=0;i<5;i++){D[i]=new double[5];}
	D[0][0]= 0;D[0][1]=10;D[0][2]=20;D[0][3]=20;D[0][4]=30;
	D[1][0]=10;D[1][1]= 0;D[1][2]=10;D[1][3]=10;D[1][4]=20;
	D[2][0]=20;D[2][1]=10;D[2][2]= 0;D[2][3]=20;D[2][4]=10;
	D[3][0]=20;D[3][1]=10;D[3][2]=20;D[3][3]= 0;D[3][4]=10;
	D[4][0]=30;D[4][1]=20;D[4][2]=10;D[4][3]=10;D[4][4]= 0;

	// random starting positions
	srand(time(0));
	Matrix2xN X;
	X.resize(5);
	for(unsigned int i=0;int(i)<X.rows();i++){
		X.Set(i,Vertex(rand()%p.maxc, rand()%p.maxc));
	}

	return SetEdges(ProblemSet(X,D));
}

#include "HundredNodeTriangulation.h"
ProblemSet Triangles100(GenParams &){
	HundredNodeTriangulation hnt;
	//var delta = new DenseVector(new double[] { 1, 1 });

	ProblemSet T;
	T.D=new double*[100];
	for(int i=0;i<100;i++){
		T.D[i]=new double[100];
		for(int j=0;j<100;j++){
			T.D[i][j]=hnt.D[i][j];
		}
	}
	T.X=hnt.X;

	T.AddLock(0,hnt.X(0,0),hnt.X(1,0));
	T.AddLock(1,hnt.X(0,1),hnt.X(1,1));
	T.AddLock(2,hnt.X(0,2),hnt.X(1,2));
	T.AddLock(3,hnt.X(0,3),hnt.X(1,3));
	//T.AddLock(50,hnt.X(0,50),hnt.X(1,50));

	return SetEdges(T);
}

using namespace TriangleProjection;
ProblemSet CTriangles40(GenParams & params){
	ProblemSet p;
	p.Init(40);
	for(int i=0;i<8;i++){
		for(int j=0;j<5;j++){
			p.X.Set(i*5+j,Vertex(i*100,j*100));
		}
	}
	p.AddLock(0,p.X(0,0),p.X(1,0));
	p.AddLock(4,p.X(0,4),p.X(1,4));
	p.AddLock(35,p.X(0,35),p.X(1,35));
	p.AddLock(39,p.X(0,39),p.X(1,39));

	for(int i=0;i<8;i++){
		for(int j=0;j<5;j++){
			if(j!=4 && i!=7){
				Triangle t((i+1)*5+j, i*5+j, i*5+j+1);
				QVector<QVector2D> tmp;t.AddToQVector(p.X,tmp);
				Q_ASSERT(IsClockwiseTriangle(tmp));
				p.D[t.a][t.b]=(p.X.GetV(t.a)-p.X.GetV(t.b)).length();
				p.D[t.b][t.c]=(p.X.GetV(t.b)-p.X.GetV(t.c)).length();
				p.D[t.c][t.a]=(p.X.GetV(t.c)-p.X.GetV(t.a)).length();
				p.D[t.b][t.a]=p.D[t.a][t.b];
				p.D[t.c][t.b]=p.D[t.b][t.c];
				p.D[t.a][t.c]=p.D[t.c][t.a];
				p.AddTriangle(t);
			}
			if(j!=0 && i!=0){
				Triangle t(i*5+j, i*5+j-1, (i-1)*5+j);
				QVector<QVector2D> tmp;t.AddToQVector(p.X,tmp);
				Q_ASSERT(IsClockwiseTriangle(tmp));
				p.D[t.a][t.b]=(p.X.GetV(t.a)-p.X.GetV(t.b)).length();
				p.D[t.b][t.c]=(p.X.GetV(t.b)-p.X.GetV(t.c)).length();
				p.D[t.c][t.a]=(p.X.GetV(t.c)-p.X.GetV(t.a)).length();
				p.D[t.b][t.a]=p.D[t.a][t.b];
				p.D[t.c][t.b]=p.D[t.b][t.c];
				p.D[t.a][t.c]=p.D[t.c][t.a];
				p.AddTriangle(t);
			}
		}
	}
	TextureLoader::FromFile(params.image,p.X);
	return SetEdges(p);
}


ProblemSet CTrianglesMxN(GenParams & params){
	ProblemSet p;
	int N=params.grid;
	int M=params.grid;
	p.Init(M*N);
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			p.X.Set(i*N+j,Vertex(i*100,j*100));
		}
	}
	p.AddLock(0,p.X(0,0),p.X(1,0));
	p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			if(j!=N-1 && i!=M-1){
				Triangle t((i+1)*N+j, i*N+j, i*N+j+1);
				QVector<QVector2D> tmp;t.AddToQVector(p.X,tmp);
				Q_ASSERT(IsClockwiseTriangle(tmp));
				p.D[t.a][t.b]=(p.X.GetV(t.a)-p.X.GetV(t.b)).length();
				p.D[t.b][t.c]=(p.X.GetV(t.b)-p.X.GetV(t.c)).length();
				p.D[t.c][t.a]=(p.X.GetV(t.c)-p.X.GetV(t.a)).length();
				p.D[t.b][t.a]=p.D[t.a][t.b];
				p.D[t.c][t.b]=p.D[t.b][t.c];
				p.D[t.a][t.c]=p.D[t.c][t.a];
				p.AddTriangle(t);
			}
			if(j!=0 && i!=0){
				Triangle t(i*N+j, i*N+j-1, (i-1)*N+j);
				QVector<QVector2D> tmp;t.AddToQVector(p.X,tmp);
				Q_ASSERT(IsClockwiseTriangle(tmp));
				p.D[t.a][t.b]=(p.X.GetV(t.a)-p.X.GetV(t.b)).length();
				p.D[t.b][t.c]=(p.X.GetV(t.b)-p.X.GetV(t.c)).length();
				p.D[t.c][t.a]=(p.X.GetV(t.c)-p.X.GetV(t.a)).length();
				p.D[t.b][t.a]=p.D[t.a][t.b];
				p.D[t.c][t.b]=p.D[t.b][t.c];
				p.D[t.a][t.c]=p.D[t.c][t.a];
				p.AddTriangle(t);
			}
		}
	}
	TextureLoader::FromFile(params.image,p.X);
	return SetEdges(p);
}
*/
ProblemSet DTrianglesMxN(GenParams & params){
	ProblemSet p;
	int N=params.grid;
	int M=params.grid;
	QVector<QVector2D> points;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D(i*(900/(M-1)),j*(900/(N-1))));
		}
	}
	p.mesh->SetVertices(points);
	//p.AddLock(0,p.X(0,0),p.X(1,0));
	//p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	//p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	//p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	TextureLoader::FromFile(params,p);
	return p;
}


/*ProblemSet DRealMxN(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	srand(params.seed);
	int N=params.grid;
	int M=params.grid;
	int V=params.count;
	int size=(M*N)+V;


	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D(i*100,j*100));
		}
	}
	//p.AddLock(0,p.X(0,0),p.X(1,0));
	//p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	//p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	//p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	//add "Real" points
	for(int i=0;i<V;i++){
		Vertex v=QVector2D(rand()%((N-1)*100),rand()%((M-1)*100));
		v.SetReal(true);
		points.append(v);
	}

	p.mesh->SetVertices(points);
	//add "Real" dissimilarities
	for(int i=M*N;i<size;i++){
		for(int j=M*N;j<size;j++){
			if(j==i){continue;}
			double dis=(p.X.GetVByIndex(i)-p.X.GetVByIndex(j)).length();
			double offset=dis/(2+rand()%4);
			if(rand()%2==0){dis-=offset;}
			else{dis+=offset;}
			p.mesh->SetDis(i,j,dis);
		}
	}

	TextureLoader::FromFile(params.image,p.X);
	return p;
}*/

#define ADDREAL(a,b,c,d) v.SetName(d);v.setX(a*1.09);v.setY(b*0.97);v.value=c;points.append(v);growth[i]=c;i++;v.RenewID();


ProblemSet SmallTest(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	srand(params.seed);
	params.grid=2;
	int M=params.grid;
	int N=params.grid;
	int V=8;
	int size=(M*N)+V;

	p.hasAssociatedValues=true;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D(i*(900/(M-1)),j*(900/(N-1))));
		}
	}
	//p.AddLock(0,p.X(0,0),p.X(1,0));
	//p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	//p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	//p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	//add "Real" points
	Vertex v;
	v.SetReal(true);
	int i=0;
	double growth[V];
	ADDREAL(681  , 52  , 1.1, "Hobart (1.1)"); //Hobart
	ADDREAL(638  ,188.2, 6.8, "Melbourne (6.8)"); //Melbourne
	ADDREAL(702  ,247  , 0.6, "Canberra (0.6)"); //Canberra
	ADDREAL(737.2,272  ,11.4, "Sydney (11.4)"); //Sydney
	ADDREAL(774  ,449.5, 4.1, "Brisbane (4.1)"); //Brisbane
	ADDREAL(109-10  ,331  , 8.6, "Perth (8.6)"); //Perth
	ADDREAL(519.2,269  , 1.0, "Adelaide (1.0)"); //Adelaide
	ADDREAL(378  ,828  , 6.0, "Darwin (6.0)"); //Darwin

	p.mesh->SetVertices(points);
	//add "Real" dissimilarities
	for(int i=M*N;i<size;i++){
		for(int j=M*N;j<size;j++){
			if(j==i){continue;}
			double dis=fabs(growth[i-M*N]-growth[j-M*N])*100;
			p.mesh->SetDis(i,j,dis);
		}
	}
	double m[6][6];
	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){
			m[i][j]=p.mesh->Dis(i,j);
		}
	}
	TextureLoader::FromFile(params,p);

	return p;
}


ProblemSet DAustralia(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	srand(params.seed);
	int N=params.grid;
	int M=params.grid;
	int V=8;
	int size=(M*N)+V;

	p.hasAssociatedValues=true;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D(i*(900/(M-1)),j*(900/(N-1))));
		}
	}
	//p.AddLock(0,p.X(0,0),p.X(1,0));
	//p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	//p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	//p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	//add "Real" points
	Vertex v;
	int idOffset=v.ID();
	v.SetReal(true);
	int i=0;
	double growth[8];
	ADDREAL(681  , 52  , 1.1, "Hobart (1.1)"); //Hobart
	ADDREAL(638  ,188.2, 6.8, "Melbourne (6.8)"); //Melbourne
	ADDREAL(702  ,247  , 0.6, "Canberra (0.6)"); //Canberra
	ADDREAL(737.2,272  ,11.4, "Sydney (11.4)"); //Sydney
	ADDREAL(774  ,449.5, 4.1, "Brisbane (4.1)"); //Brisbane
	ADDREAL(109-10  ,331  , 8.6, "Perth (8.6)"); //Perth
	ADDREAL(519.2,269  , 1.0, "Adelaide (1.0)"); //Adelaide
	ADDREAL(378  ,828  , 6.0, "Darwin (6.0)"); //Darwin

	p.mesh->SetVertices(points);
	//add "Real" dissimilarities
	for(int i=0;i<V;i++){
		for(int j=0;j<V;j++){
			if(j==i){continue;}
			double dis=fabs(growth[i]-growth[j])*100;
			p.mesh->SetDis(i+idOffset,j+idOffset,dis);
		}
	}

	TextureLoader::FromFile(params,p);
	return p;
}

ProblemSet DAustraliaM(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	srand(params.seed);
	int N=params.grid;
	int M=params.grid;
	int V=8;
	int size=(M*N)+V;

	p.hasAssociatedValues=true;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D(i*(900/(M-1)),j*(900/(N-1))));
		}
	}
	//p.AddLock(0,p.X(0,0),p.X(1,0));
	//p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	//p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	//p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	p.mesh->SetVertices(points);
	points.clear();
	//add "Real" points
	ProblemSet pReal;
	Vertex v;
	int idOffset=v.ID();
	v.SetReal(true);
	int i=0;
	double growth[8];
	ADDREAL(681  , 52  , 1.1, "Hobart (1.1)"); //Hobart
	ADDREAL(638  ,188.2, 6.8, "Melbourne (6.8)"); //Melbourne
	ADDREAL(702  ,247  , 0.6, "Canberra (0.6)"); //Canberra
	ADDREAL(737.2,272  ,11.4, "Sydney (11.4)"); //Sydney
	ADDREAL(774  ,449.5, 4.1, "Brisbane (4.1)"); //Brisbane
	ADDREAL(109-10  ,331  , 8.6, "Perth (8.6)"); //Perth
	ADDREAL(519.2,269  , 1.0, "Adelaide (1.0)"); //Adelaide
	ADDREAL(378  ,828  , 6.0, "Darwin (6.0)"); //Darwin

	pReal.mesh->SetVertices(points);

	//add "Real" dissimilarities
	for(int i=0;i<V;i++){
		for(int j=0;j<V;j++){
			if(j==i){continue;}
			double dis=fabs(growth[i]-growth[j])*100;
			pReal.mesh->SetDis(i+idOffset,j+idOffset,dis);
		}
	}

	AddSolutionAsMovement(pReal,p);

	TextureLoader::FromFile(params,p);
	return p;
}

#define ADDREALP(a,b,c) v.SetName(c);v.setX((a*1.1)/9*(params.grid-1));v.setY((b*0.97)/9*(params.grid-1));points.append(v);v.RenewID();
//#define ADDREALP(a,b) v.setX(a);v.setY(b);points.append(v);v.RenewID();
#define DIS(a,b,c) p.mesh->SetDis(idOffset+a,idOffset+b,c);

ProblemSet DRoads(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	srand(params.seed);
	int N=params.grid;
	int M=params.grid;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D(i*(900/(M-1)),j*(900/(N-1))));
		}
	}
	//p.AddLock(0,p.X(0,0),p.X(1,0));
	//p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	//p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	//p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	//add "Real" points
	Vertex v;
	int idOffset=v.ID();
	v.SetReal(true);
	ADDREALP(681, 52, "Hobart"); //Hobart
	ADDREALP(638,188.2, "Melbourne"); //Melbourne
	ADDREALP(702,247, "Canberra"); //Canberra
	ADDREALP(737.2,272, "Sydney"); //Sydney
	ADDREALP(774,449.5, "Brisbane"); //Brisbane
	ADDREALP(109-10,331, "Perth"); //Perth
	ADDREALP(519.2,269, "Adelaide"); //Adelaide
	ADDREALP(378,828, "Darwin"); //Darwin

	p.mesh->SetVertices(points);
	//add "Real" dissimilarities


	DIS(0,1,23*60+39); //733
	DIS(0,2,30*60); //1396
	DIS(0,3,32*60); //1611
	DIS(0,4,42*60); //2419
	DIS(0,5,59*60); //4143
	DIS(0,6,31*60); //1452
	DIS(0,7,65*60); //4478
	DIS(1,2,6*60+34); //663
	DIS(1,3,8*60+38); //878
	DIS(1,4,18*60+12); //1686
	DIS(1,5,35*60); //3418
	DIS(1,6,7*60+49); //727
	DIS(1,7,42*60); //3753
	DIS(2,3,2*60+56); //286
	DIS(2,4,12*60+36); //1200
	DIS(2,5,38*60); //3719
	DIS(2,6,12*60+8); //1160
	DIS(2,7,45*60); //3940
	DIS(3,4,10*60+3); //924
	DIS(3,5,40*60); //3935
	DIS(3,6,14*60+10); //1376
	DIS(3,7,45*60); //3929
	DIS(4,5,45*60); //4344
	DIS(4,6,21*60+39); //2024
	DIS(4,7,39*60); //3425
	DIS(5,6,28*60); //2693
	DIS(5,7,43*60); //4032
	DIS(6,7,34*60); //3029


	TextureLoader::FromFile(params,p);
	return p;
}

ProblemSet DMove(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	srand(params.seed);
	int N=params.grid;
	int M=params.grid;
	int V=params.count;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D(i*(900/(M-1)),j*(900/(N-1))));
		}
	}
	//p.AddLock(0,p.X(0,0),p.X(1,0));
	//p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
	//p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
	//p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

	//add "moving" points
	p.mesh->SetVertices(points);

	srand(params.seed);
	for(int i=0;i<V;i++){
		QVector2D end(100*(double((M-1))/V*i+0.5),(N-1)*0.5*100);
		Vertex * pos= new Vertex(end + QVector2D((rand()%100-50)*(M/10),(rand()%200-100))*N/10);
		pos->SetReal(true);
		p.mesh->AddVertex(pos);
		p.moves.append(MovementData(pos,*pos,end));
	}

	p.mesh->Triangulate();

	TextureLoader::FromFile(params,p);
	return p;
}

#include "londonbikes.cpp"
ProblemSet LondonBikes(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	int M=params.grid;
	int N=params.grid;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D((i-0.5)*(15600/(M-2)),(j-0.5)*(7008/(N-2))));
		}
	}
	p.mesh->SetVertices(points);
	LondenBikes(p);
	p.mesh->Triangulate();

	TextureLoader::FromFile(params,p);
	return p;
}


#define ADDREALM(a,b,c,d) v=new Vertex(a,b,true);p.mesh->AddVertex(v);p.moves.append(MovementData(v,(*v),QVector2D(c,d)));
ProblemSet LondonRiver(GenParams & params){
	ProblemSet p;
	QVector<Vertex> points;
	int M=params.grid;
	int N=params.grid;

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			points.append(QVector2D((i-0.5)*(15600/(M-2)),(j-0.5)*(7008/(N-2))));
		}
	}

	p.mesh->SetVertices(points);

	Vertex *v;
	ADDREALM(3778.4, -164.685,3201.4, 72.9031);
	ADDREALM(4932.4, 101.187,3241, 1300.44);
	ADDREALM(6374.9, 191.697,3218.38, 2182.91);
	ADDREALM(6980.18, 1209.93,4513.79, 2182.91);
	ADDREALM(7166.86, 2607.17,6929.27, 2194.22);
	ADDREALM(8219.03, 3206.8,8615.01, 2222.51);
	ADDREALM(10538.3, 2805.16,10577.9, 2194.22);
	ADDREALM(12546.5, 2906.99,12546.5, 2194.22);
	ADDREALM(13989, 2748.6,14424.6, 2222.51);
	ADDREALM(14435.9, 723.441,14571.7, 1594.6);
	ADDREALM(16251.8, 1328.72,14652.6, -40.335);


	TextureLoader::FromFile(params,p);
	return p;
}

ProblemSet DNuts(GenParams & params, std::function<void(QVector<QVector2D> &reals, int,double,double, double,QVector<Vertex> &)> AddHelperPoints){
    QFile areaData("../data/nuts1/NUTS1_AreaData.csv");
    if (!areaData.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "File error";
    }

    areaData.readLine(); // skip header
    QMap<QString, std::tuple<double, double>> areaCentroids;
    typedef std::numeric_limits<double> limitsd;
    double minX = limitsd::max(),
           maxX = limitsd::lowest(),
           minY = limitsd::max(),
           maxY = limitsd::lowest();
    while (!areaData.atEnd()) {
        QString line = areaData.readLine();
        QStringList cells = line.split(',');
        double x = cells[3].toDouble() * 100;
        double y = cells[4].toDouble() * 100;
        areaCentroids[cells.first()] = std::make_tuple(x, y);
        minX = min(minX, x);
        minY = min(minY, y);
        maxX = max(maxX, x);
        maxY = max(maxY, y);
    }

    QFile file("../data/nuts1/OACNUTS1_RNG.csv");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "File error";
    }

    QVector<QVector<double>> data;
    int fields;
    file.readLine(); // skip header
    while (!file.atEnd()) {
        QString line = file.readLine();
        QStringList cells = line.split(',');
        cells.removeFirst();
        fields = cells.size();
        QVector<double> row;
        foreach (QString cell, cells) {
            row.append(cell.toDouble());
        }
        data.append(row);
    }

    ProblemSet p;
    QVector<Vertex> points;
    srand(params.seed);

		double aspectRatio=13.0/7.0;
		int V=areaCentroids.size();
		//int size=(M*N)+V;
    double width = maxX - minX;
    double height = maxY - minY;


		p.hasAssociatedValues=true;
				//p.AddLock(0,p.X(0,0),p.X(1,0));
    //p.AddLock(N-1,p.X(0,N-1),p.X(1,N-1));
    //p.AddLock(N*(M-1),p.X(0,N*(M-1)),p.X(1,N*(M-1)));
    //p.AddLock(M*N-1,p.X(0,M*N-1),p.X(1,M*N-1));

    //add "Real" points

    auto addrealp=[&](double a, double b, QString c) {
        Vertex v;
				//idOffset = v.ID(); //I used this before to store the first ID which identified a real data point so I could set their dissimilarities manually, as it is now it stores the id of the last real vertex which was added
        v.SetReal(true);
        v.SetName(c);
        v.setX(170 + 0.425 * (a-minX));
        v.setY(70 + 0.645 * (b-minY));
        points.append(v);
				//v.RenewID();//only needed when you reuse a vertex (as in, the actual variable), the id gets incremented automatically in the constructor, so now you are incrementing it twice per datapoint
    };

		QVector<QVector2D> reals;
    foreach (QString areaName, areaCentroids.keys()) {
        double x, y;
        std::tie<double,double>(x,y) = areaCentroids.value(areaName);
				addrealp(x, y, areaName);
				reals.append(points.last());
    }

		AddHelperPoints(reals,params.grid,aspectRatio,width,height,points);

		p.mesh->SetVertices(points);
		int idOffset=points.first().ID();
    //add "Real" dissimilarities
    for(int i=0;i<V;i++){
        for(int j=0;j<V;j++){
            if(j==i){continue;}
            double dis = 0;
            for(int k = 0; k < fields; ++k) {
                double dij = data[i][k] - data[j][k];
                dis += dij*dij;
            }
            dis = 100 * sqrt(dis);
            p.mesh->SetDis(i+idOffset,j+idOffset,dis);
        }
    }

    TextureLoader::FromFile(params,p);
    return p;
}

ProblemSet DNutsGrid(GenParams &params){
	return DNuts(params,[](QVector<QVector2D> &reals, int pts, double aspectRatio,double width, double height, QVector<Vertex> &points){
		int N=aspectRatio * pts;
		int M=pts;
		for(int i=0;i<M;i++){
				for(int j=0;j<N;j++){
						points.append(QVector2D(i*(width/(M-1)),aspectRatio * j*(height/(N-1))));
				}
		}
	});
}


ProblemSet DNutsNonUniform(GenParams &params){
	return DNuts(params,[&params](QVector<QVector2D> &reals, int pts, double aspectRatio,double width, double height, QVector<Vertex> &points){
		QRectF r(QPointF(0,0),QPointF(width,height*aspectRatio));
		MeshPoints(reals,r,points,height/params.grid);
	});
}

Generator::type Generator::generators[] = {
		{"NUTS Grid", DNutsGrid},
	{"NUTS QuadTree", DNutsNonUniform},
	{"DTriangleMxN", DTrianglesMxN},
	//{"RealTest",SmallTest},
	//{"Real", DRealMxN},
	{"House prices", DAustralia},
	{"House prices - M", DAustraliaM},
	{"Travel time", DRoads},
	{"DMove", DMove},
	{"LondonBikes", LondonBikes},
	//{"London", LondonRiver},
    {"File", DFile},

	{NULL,NULL}, // Last element of list.
};

