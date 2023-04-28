#include <QDebug>
#include <QTimer>
#include <QLabel>
#include <QSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QColorDialog>
#include <QTextEdit>
#include <QComboBox>

#ifdef KDE4_FOUND
#include <ktexteditor/configinterface.h>
#include <ktexteditor/editor.h>
#include <ktexteditor/editorchooser.h>
#include <ktexteditor/view.h>
#endif

#include "config.h"

Config::Config(QWidget *parent, Qt::WindowFlags f)
	: QWidget(parent, f), vbox(new QVBoxLayout(this)) {
}

void Config::done()
{
	vbox->addStretch(1);
}

Vector3DItem::Vector3DItem(Config *parent, QString name)
	: QObject(parent), changing(false) {
	parent->vbox->addWidget(new QLabel(name+":", parent));
	QHBoxLayout *box = new QHBoxLayout();
	parent->vbox->addLayout(box);
	x = new QDoubleSpinBox(parent);
	x->setRange(-1000, 1000);
	connect(x, SIGNAL(valueChanged(double)), SLOT(setX(double)));
	box->addWidget(x);
	y = new QDoubleSpinBox(parent);
	y->setRange(-1000, 1000);
	connect(y, SIGNAL(valueChanged(double)), SLOT(setY(double)));
	box->addWidget(y);
	z = new QDoubleSpinBox(parent);
	z->setRange(-1000, 1000);
	connect(z, SIGNAL(valueChanged(double)), SLOT(setZ(double)));
	box->addWidget(z);
}

void Vector3DItem::setVector(QVector3D v) {
	value = v;
	changing = true;
	x->setValue(v.x());
	y->setValue(v.y());
	z->setValue(v.z());
	changing = false;
}

void Vector3DItem::setX(double x) {
	if (changing) return;
	value.setX(x);
	emit vectorChanged(value);
}

void Vector3DItem::setY(double y) {
	if (changing) return;
	value.setY(y);
	emit vectorChanged(value);
}

void Vector3DItem::setZ(double z) {
	if (changing) return;
	value.setZ(z);
	emit vectorChanged(value);
}

void Vector3DItem::setRange(double min, double max)
{
	x->setRange(min, max);
	y->setRange(min, max);
	z->setRange(min, max);
}

void Vector3DItem::setSingleStep(double v)
{
	x->setSingleStep(v);
	y->setSingleStep(v);
	z->setSingleStep(v);
}

RangeItem::RangeItem(Config* parent, QString name)
	: QObject(parent), changing(false) {
	parent->vbox->addWidget(new QLabel(name+":", parent));
	QHBoxLayout *box = new QHBoxLayout();
	parent->vbox->addLayout(box);
	smin = new QDoubleSpinBox(parent);
	connect(smin, SIGNAL(valueChanged(double)), SLOT(setMin(double)));
	box->addWidget(smin);
	smax = new QDoubleSpinBox(parent);
	connect(smax, SIGNAL(valueChanged(double)), SLOT(setMax(double)));
	box->addWidget(smax);
}

void RangeItem::setValue(double vmin, double vmax)
{
	cmin = vmin;
	cmax = vmax;
	changing = true;
	smin->setValue(vmin);
	smax->setValue(vmax);
	changing = false;
}

void RangeItem::setMin(double v)
{
	if (changing) return;
	cmin = v;
	emit valueChanged(cmin,cmax);
}

void RangeItem::setMax(double v)
{
	if (changing) return;
	cmax = v;
	emit valueChanged(cmin,cmax);
}

void RangeItem::setRange(double min, double max) {
	smin->setRange(min, max);
	smax->setRange(min, max);
}

void RangeItem::setSingleStep(double v) {
	smin->setSingleStep(v);
	smax->setSingleStep(v);
}

void RangeItem::setDecimals(int v) {
	smin->setDecimals(v);
	smax->setDecimals(v);
}

DoubleItem::DoubleItem(Config* parent, QString name, double * var)
    : QObject(parent), changing(false), var_ptr(var) {
    QHBoxLayout *box = new QHBoxLayout();
    parent->vbox->addLayout(box);
    box->addWidget(new QLabel(name+":", parent));
    spinner = new QDoubleSpinBox(parent);
    spinner->setRange(-1000, 1000);
    if (var_ptr) spinner->setValue(*var_ptr);
    connect(spinner, SIGNAL(valueChanged(double)), SLOT(setValueProxy(double)));
    box->addWidget(spinner);
}

void DoubleItem::setValue(double v) {
    changing = true;
    spinner->setValue(v);
    changing = false;
}

void DoubleItem::setValueProxy(double v) {
    if (var_ptr) *var_ptr=v;
    if (changing) return;
    emit valueChanged(v);
}

void DoubleItem::setRange(double min, double max) {
    spinner->setRange(min, max);
}

void DoubleItem::setSingleStep(double v) {
    spinner->setSingleStep(v);
}

void DoubleItem::setDecimals(int v) {
    spinner->setDecimals(v);
}

IntegerItem::IntegerItem(Config* parent, QString name, int * var)
    : QObject(parent), changing(false), var_ptr(var) {
    QHBoxLayout *box = new QHBoxLayout();
    parent->vbox->addLayout(box);
    box->addWidget(new QLabel(name+":", parent));
    spinner = new QSpinBox(parent);
    spinner->setRange(-1000, 1000);
    if (var_ptr) spinner->setValue(*var_ptr);
    connect(spinner, SIGNAL(valueChanged(int)), SLOT(setValueProxy(int)));
    box->addWidget(spinner);
}

void IntegerItem::setValue(int v) {
    changing = true;
    spinner->setValue(v);
    changing = false;
}

void IntegerItem::setValueProxy(int v) {
    if (var_ptr) *var_ptr=v;
    if (changing) return;
    emit valueChanged(v);
}

void IntegerItem::setRange(int min, int max) {
    spinner->setRange(min, max);
}

void IntegerItem::setSingleStep(int v) {
    spinner->setSingleStep(v);
}

OutputItem::OutputItem(Config * parent): QObject(parent){
	QHBoxLayout *box = new QHBoxLayout();
	parent->vbox->addLayout(box);
	txt = new QTextEdit(parent);
	txt->setReadOnly(true);
	txt->setFontFamily("Monospace");
	box->addWidget(txt);
}

void OutputItem::appendText(QString v){
	txt->append(v);
	txt->moveCursor(QTextCursor::End);
}

LabelItem::LabelItem(Config* parent)
	: QObject(parent) {
	QHBoxLayout *box = new QHBoxLayout();
	parent->vbox->addLayout(box);
	label = new QLabel(parent);
	box->addWidget(label);
}

void LabelItem::setValue(QString v) {
	label->setText(v);
}

void LabelItem::appendText(QString v){
	QString txt = label->text().append(v);
	label->setText(txt);
}

void LabelItem::setValue(double v)
{
	QString s;
	s.sprintf("%.2lf", v);
	label->setText(s);
}

void LabelItem::setValue(int v)
{
	QString s;
	s.sprintf("%d", v);
	label->setText(s);
}

NamedLabelItem::NamedLabelItem(Config* parent, QString name)
	: QObject(parent) {
	QHBoxLayout *box = new QHBoxLayout();
	parent->vbox->addLayout(box);
	box->addWidget(new QLabel(name+":", parent));
	label = new QLabel(parent);
	box->addWidget(label);
}

void NamedLabelItem::setValue(QString v) {
	label->setText(v);
}

void NamedLabelItem::setValue(double v)
{
	QString s;
	s.sprintf("%.2lf", v);
	label->setText(s);
}

void NamedLabelItem::setValue(int v)
{
	QString s;
	s.sprintf("%d", v);
	label->setText(s);
}

BoolItem::BoolItem(Config* parent, QString name, bool * var)
	: QObject(parent), changing(false), var_ptr(var) {
	cbox = new QCheckBox(name, parent);
    if (var_ptr) cbox->setChecked(*var_ptr);
	connect(cbox, SIGNAL(toggled(bool)), SLOT(setValueProxy(bool)));
	setValue(*var);
	parent->vbox->addWidget(cbox);
}

void BoolItem::setValue(bool v)
{
	changing = true;
	cbox->setChecked(v);
	changing = false;
}

void BoolItem::setValueProxy(bool v)
{
    if (var_ptr) *var_ptr = v;
	if (changing) return;
	emit valueChanged(v);
}

ColorItem::ColorItem(Config* parent, QString name, bool alpha)
	: QObject(parent), alpha(alpha) {
	QHBoxLayout *box = new QHBoxLayout();
	parent->vbox->addLayout(box);
	box->addWidget(new QLabel(name+":", parent));
	button = new QPushButton(parent);
	col = new QColorDialog(parent);
	connect(button, SIGNAL(clicked(bool)), col, SLOT(open()));
	connect(col, SIGNAL(colorSelected(QColor)), SIGNAL(valueChanged(QColor)));
	connect(col, SIGNAL(colorSelected(QColor)), SLOT(setValueProxy(QColor)));
	if (alpha) col->setOption(QColorDialog::ShowAlphaChannel);
	box->addWidget(button);
}

void ColorItem::setValue(QColor v)
{
	col->setCurrentColor(v);
	setValueProxy(v);
}

void ColorItem::setValueProxy(QColor v)
{
	button->setStyleSheet("* { background-color: "+v.name()+"; color: "+(v.valueF()>.5?"black":"white")+" }");
	if (alpha) {
		QString s;
		s.sprintf("%.0f%%", v.alphaF()*100);
		button->setText(s);
	}
}

FileChooserItem::FileChooserItem(Config* parent, QString name, QString *filename, QString filter, QString defaultDir)
	: QObject(parent), filename(filename), name(name), filter(filter), defaultDir(defaultDir){
	btn = new QPushButton(parent);
	btn->setText(name);
	if(*filename!=""){
		QFileInfo info=(*filename);
		btn->setText(name + ": " + info.fileName());
	}
	parent->vbox->addWidget(btn);
	connect(btn, SIGNAL(clicked()), SLOT(OpenFile()));
}

TextBoxItem::TextBoxItem(Config* parent, QString name)
	: QObject(parent), changing(false) {
	parent->vbox->addWidget(new QLabel(name+":", parent));
#ifdef KDE4_FOUND
	edit = KTextEditor::EditorChooser::editor()->createDocument(parent);
	view = edit->createView(parent);
	parent->vbox->addWidget(view);
	connect(edit, SIGNAL(textChanged(KTextEditor::Document*)), SLOT(changed()));
	KTextEditor::ConfigInterface *iface = qobject_cast<KTextEditor::ConfigInterface*>( view );
	if (iface) {
		def_color = iface->configValue("background-color");
	}
#else
	edit = new QTextEdit(parent);
	connect(edit, SIGNAL(textChanged()), SLOT(changed()));
	parent->vbox->addWidget(edit);
	edit->setAcceptRichText(false);
#endif
	timer = new QTimer(parent);
	timer->setSingleShot(true);
	connect(timer, SIGNAL(timeout()), SLOT(send()));
}

void TextBoxItem::setLanguage(QString /*lang*/)
{
#ifdef KDE4_FOUND
	edit->setHighlightingMode(lang);
#endif
}

void TextBoxItem::changed()
{
	if (changing) return;
#ifdef KDE4_FOUND
	KTextEditor::ConfigInterface *iface = qobject_cast<KTextEditor::ConfigInterface*>( view );
	if (iface) {
		iface->setConfigValue("background-color", QColor("#ffffcc"));
	}
#else
	edit->setStyleSheet("* { background-color: #ffffcc; }");
#endif
	timer->start(500);
}

void TextBoxItem::send()
{
#ifdef KDE4_FOUND
	KTextEditor::ConfigInterface *iface = qobject_cast<KTextEditor::ConfigInterface*>( view );
	if (iface) {
		iface->setConfigValue("background-color", def_color);
	}
	emit valueChanged(edit->text());
#else
	edit->setStyleSheet("");
	emit valueChanged(edit->toPlainText());
#endif
}

void TextBoxItem::setValue(QString v)
{
	changing = true;
	edit->setText(v);
	changing = false;
}

OptionItem::OptionItem(Config* parent, QString name, int * var)
    : QObject(), changing(false), var_ptr(var)
{
	QHBoxLayout *box = new QHBoxLayout();
	parent->vbox->addLayout(box);
	box->addWidget(new QLabel(name+":", parent));
	list = new QComboBox(parent);
	box->addWidget(list);
	connect(list, SIGNAL(currentIndexChanged(int)), SLOT(setValueProxy(int)));
}

int OptionItem::addItem(QString text)
{
    changing=true;
	list->addItem(text);
    changing=false;
    if (var_ptr && *var_ptr==list->count()-1) list->setCurrentIndex(*var_ptr);
	return list->count()-1;
}

int OptionItem::addSeparator()
{
	list->insertSeparator(list->count());
	return list->count()-1;
}

void OptionItem::setToOther()
{
	changing=true;
	list->setCurrentIndex(other);
	changing=false;
}

void OptionItem::setValue(int i)
{
	changing=true;
	list->setCurrentIndex(i);
	changing=false;
}

void OptionItem::setValueProxy(int i)
{
	if (changing) return;
    if (var_ptr) *var_ptr = i;
	emit valueChanged(i);
}

ButtonItem::ButtonItem(Config *parent, QString caption): QPushButton(caption, parent) {
    parent->vbox->addWidget(this);
}

void ButtonItem::disable() {
    setEnabled(false);
}

void ButtonItem::enable() {
    setEnabled(true);
}
