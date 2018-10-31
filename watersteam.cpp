#include "watersteam.h"
#include "ui_watersteam.h"
#include "qmath.h"
#include <QList>
#include <QDebug>

double t_Conversion_K_to_Any(int i, double j)
{
    if (i==0) {
        j = j-273.15;
    } else if (i==1) {
        j = (j-273.15)*(9.00/5.00)+32;
    } else if (i==2) {
        j = j;
    } else {
        j = j*(9.00/5.00);
    }
    return(j);
}

double t_Conversion_Any_to_K(int i, double j)
{
    if (i==0) {
        j = j+273.15;
    } else if (i==1) {
        j = (j-32)*(5.00/9.00)+273.15;
    } else if (i==2) {
        j = j;
    } else {
        j = j*(5.00/9.00);
    }
    return(j);
}

WaterSteam::WaterSteam(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::WaterSteam)
{
    ui->setupUi(this);
}

WaterSteam::~WaterSteam()
{
    delete ui;
}
