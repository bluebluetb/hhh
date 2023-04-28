#ifndef CONFIG_H
#define CONFIG_H

#include <QtCore/QObject>
#include <QtCore/QVariant>
#include <QWidget>
#include <QVBoxLayout>
#include <QVector3D>
#include <QTextEdit>
#include <QPushButton>
#include <QFileDialog>
#include <QLabel>

class QSpinBox;
namespace KTextEditor {
class Document;
class View;
}

class QComboBox;
class QTextEdit;
class QColorDialog;
class QPushButton;
class QCheckBox;
class QLabel;
class QDoubleSpinBox;
class Config;

class Config : public QWidget {
    Q_OBJECT;

public:
    explicit Config(QWidget *parent = 0, Qt::WindowFlags f = 0);
    virtual void done();

    QVBoxLayout *vbox;
};

class Vector3DItem : public QObject {
	Q_OBJECT;
private:
	QVector3D value;
	QDoubleSpinBox *x, *y, *z;
	bool changing; // Used to prevent a value changed event when setting the value.
public:
	explicit Vector3DItem(Config * parent, QString name);
public slots:
	virtual void setVector(QVector3D v);
	virtual void setRange(double min, double max);
	virtual void setSingleStep(double v);
private slots:
	virtual void setX(double x);
	virtual void setY(double y);
	virtual void setZ(double z);
signals:
	void vectorChanged(QVector3D);
};

class RangeItem : public QObject {
	Q_OBJECT;
private:
	double cmin, cmax;
	QDoubleSpinBox *smin;
	QDoubleSpinBox *smax;
	bool changing; // Used to prevent a value changed event when setting the value.
public:
	explicit RangeItem(Config * parent, QString name);
public slots:
	virtual void setRange(double min, double max);
	virtual void setSingleStep(double v);
	virtual void setDecimals(int v);
	virtual void setValue(double vmin, double vmax);
private slots:
	virtual void setMin(double v);
	virtual void setMax(double v);
signals:
	void valueChanged(double,double);
};

class DoubleItem : public QObject {
	Q_OBJECT;
private:
	QDoubleSpinBox *spinner;
	bool changing; // Used to prevent a value changed event when setting the value.
	double * var_ptr;
public:
	explicit DoubleItem(Config * parent, QString name, double * var=NULL);
public slots:
	virtual void setRange(double min, double max);
	virtual void setSingleStep(double v);
	virtual void setDecimals(int v);
	virtual void setValue(double v);
private slots:
	virtual void setValueProxy(double v);
signals:
	void valueChanged(double);
};

class IntegerItem : public QObject {
    Q_OBJECT;
private:
    QSpinBox *spinner;
    bool changing; // Used to prevent a value changed event when setting the value.
    int * var_ptr;
public:
    explicit IntegerItem(Config * parent, QString name, int * var=NULL);
public slots:
    virtual void setRange(int min, int max);
    virtual void setSingleStep(int v);
    virtual void setValue(int v);
private slots:
    virtual void setValueProxy(int v);
signals:
    void valueChanged(int);
};

class OutputItem : public QObject {
	Q_OBJECT;
private:
	QTextEdit * txt;
public:
	explicit OutputItem(Config *parent);
public slots:
	virtual void appendText(QString v);
};

class LabelItem : public QObject {
	Q_OBJECT;
private:
	QLabel * label;
public:
	explicit LabelItem(Config *parent);
public slots:
	virtual void setValue(QString v);
	virtual void appendText(QString v);
	virtual void setValue(double v);
	virtual void setValue(int v);
	QLabel * getLabel(){return label;}
};

class NamedLabelItem : public QObject {
	Q_OBJECT;
private:
	QLabel * label;
public:
	explicit NamedLabelItem(Config *parent, QString name);
public slots:
	virtual void setValue(QString v);
	virtual void setValue(double v);
	virtual void setValue(int v);
};

class BoolItem : public QObject {
	Q_OBJECT;
private:
	QCheckBox * cbox;
	bool changing; // Used to prevent a value changed event when setting the value.
    bool * var_ptr;
public:
	explicit BoolItem(Config * parent, QString name, bool * var=NULL);
public slots:
	virtual void setValue(bool v);
private slots:
	virtual void setValueProxy(bool v);
signals:
	void valueChanged(bool);
};

class ColorItem : public QObject {
	Q_OBJECT;
private:
	QPushButton * button;
	QColorDialog * col;
	bool alpha;
public:
	explicit ColorItem(Config * parent, QString name, bool alpha=false);
public slots:
	virtual void setValue(QColor v);
private slots:
	virtual void setValueProxy(QColor v);
signals:
	void valueChanged(QColor);
};

class FileChooserItem : public QObject {
		Q_OBJECT;
	private:
		QPushButton * btn;
		QString * filename;
		QString name;
		QString filter;
		QString defaultDir;
		void UpdateButtonText(){
			QFileInfo info=(*filename);
			btn->setText(name + ": " +info.fileName());
		}
	protected slots:
		void OpenFile(){
			(*filename) = QFileDialog::getOpenFileName(0,tr("Open Image"), defaultDir, filter);
			UpdateButtonText();
		}
	public:
		void SetFile(QString file){
			(*filename)=file;
			UpdateButtonText();
		}
		explicit FileChooserItem(Config * parent, QString name, QString * filename, QString filter="All Files (*)", QString defaultDir="");

};

class TextBoxItem : public QObject {
	Q_OBJECT;
private:
#ifdef KDE4_FOUND
	KTextEditor::Document * edit;
	KTextEditor::View * view;
	QVariant def_color;
#else
	QTextEdit * edit;
#endif
	QTimer * timer;
	bool changing; // Used to prevent a value changed event when setting the value.
public:
	explicit TextBoxItem(Config * parent, QString name);
	void setLanguage(QString lang);
public slots:
	virtual void setValue(QString v);
private slots:
	void changed();
	void send();
signals:
	void valueChanged(QString);
};

class OptionItem : public QObject {
	Q_OBJECT;
private:
	QComboBox * list;
	bool changing; // Used to prevent a value changed event when setting the value.
    int * var_ptr;
public:
	explicit OptionItem(Config * parent, QString name, int * var = NULL);
	int addItem(QString text);
	int addSeparator();
	int other;
public slots:
	virtual void setValue(int i);
	virtual void setToOther();
private slots:
	virtual void setValueProxy(int i);
signals:
	void valueChanged(int);
};

class ButtonItem : public QPushButton {
    Q_OBJECT;
public:
    explicit ButtonItem(Config * parent, QString caption);
public slots:
    void enable();
    void disable();
};

#endif // CONFIG_H
