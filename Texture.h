#ifndef TEXTURE_H
#define TEXTURE_H

#include "structures.h"
#include <GL/gl.h>
#include <QImage>
#include <QGLWidget>
#include <QRgb>

class TextureLoader{
	public:
		static void FromFile(GenParams &gp, ProblemSet &p){
			QString filename=gp.image;
			QImage img;
			if(!img.load(filename)){
				DEBUG << "Error loading image - using cyan placeholder";
				img = QImage(128,128,QImage::Format_RGB32);
				img.fill(Qt::cyan);
			}
			QImage glImg = QGLWidget::convertToGLFormat(img);
			p.imgWidth=glImg.width();
			p.imgHeight=glImg.height();
			if(gp.overlayGrid){
				int deltax=p.imgWidth/(gp.overlayGridX+1);
				int deltay=p.imgHeight/(gp.overlayGridY+1);
				for(int x=1;x<=gp.overlayGridX;x++){
					for(int y=0;y<p.imgHeight;y++){
						for(int i=-1;i<=1;i++){
							QRgb color=glImg.pixel(x*deltax+i,y);
							color=qRgb(min(255,qRed(color)*1.2+50),min(255,qGreen(color)*1.2+50),min(255,qBlue(color)*1.2+50));
							glImg.setPixel(x*deltax+i,y,color);
						}
					}
				}
				for(int y=1;y<=gp.overlayGridY;y++){
					for(int x=0;x<p.imgWidth;x++){
						for(int i=-1;i<=1;i++){
							QRgb color=glImg.pixel(x,y*deltay+i);
							color=qRgb(min(255,qRed(color)+50),min(255,qGreen(color)+50),min(255,qBlue(color)+50));
							glImg.setPixel(x,y*deltay+i,color);
						}
					}
				}
			}
			GLuint tmp[1];
			glGenTextures(1, &tmp[0]);
			p.glTextureId=tmp[0];
			glBindTexture(GL_TEXTURE_2D, p.glTextureId);
			glTexImage2D( GL_TEXTURE_2D, 0, 3, glImg.width(), glImg.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, glImg.bits() );
			glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
			glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			SetCoordinates(p);
			p.hasTexture=true;
		}
	private:
		static void SetCoordinates(ProblemSet &p){
			VertexIterator it=p.mesh->VBegin();
			QVector2D min=it.GetV(),max=it.GetV();
			for (; it!=p.mesh->VEnd(); it.Next()) {
				QVector2D x=it.GetV();
				if(x.x()<min.x()){min=QVector2D(x.x(),min.y());}
				if(x.y()<min.y()){min=QVector2D(min.x(),x.y());}
				if(x.x()>max.x()){max=QVector2D(x.x(),max.y());}
				if(x.y()>max.y()){max=QVector2D(max.x(),x.y());}
			}
			p.minOrigCoord=min;
			p.maxOrigCoord=max;
			double w=max.x()-min.x();
			double h=max.y()-min.y();
			it=p.mesh->VBegin();
			for (; it!=p.mesh->VEnd(); it.Next()) {
				QVector2D x=it.GetV();
				x-=min;
				x=QVector2D(x.x()/w, x.y()/h);
				it.Vertex().tx=x.x();
				it.Vertex().ty=x.y();
			}
		}
};

/*
ProblemSet SetTestTexture(ProblemSet p){
	unsigned size=pow(2,7);
	int blocks=pow(2,3);
	int b=0,w=0;
	GLubyte * raw=new GLubyte[size * size * 3];
	for(unsigned i=0;i<size;i++){
		for(unsigned j=0;j<size;j++){
			raw[i]=0x00;
			bool black=(int(i/blocks)%2==0); //even rows
			if((int(j/blocks)%2==1)){black=!black;} //toggle if uneven column
			if(black){
				//raw[i]=0xff;
				b++;
			}else{
				raw[i*size*3+j*3]=0xff;
				raw[i*size*3+j*3+1]=0xff;
				raw[i*size*3+j*3+2]=0xff;

				w++;
			}
		}
	}
	DEBUG << "b:" <<b << "    w:"<< w;
	p.texture=new TextureData(p.X, raw,size,size);
	//delete raw;
	return p;
}*/


#endif // TEXTURE_H
