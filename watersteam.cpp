#include "watersteam.h"
#include "ui_watersteam.h"

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
