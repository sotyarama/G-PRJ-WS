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
    init();
}

WaterSteam::~WaterSteam()
{
    delete ui;
}

//[REGION - 4]
void WaterSteam::reg_4_saturation_Temperature()
{
    double t_Ref = t_Reference_region4;
    double p_Ref = p_Reference_region4;
    double pressure = p_Input;
    double beta = pow((pressure/p_Ref),0.25);

    double E = pow(beta,2)+n_region4[3]*beta+n_region4[6];
    double F = n_region4[1]*pow(beta,2)+n_region4[4]*beta+n_region4[7];
    double G = n_region4[2]*pow(beta,2)+n_region4[5]*beta+n_region4[8];
    double D = (2.00*G)/(-F-pow((pow(F,2)-4.00*E*G),0.50));

    double temperature = t_Ref*(n_region4[10]+D-pow((pow((n_region4[10]+D),2)-4.00*(n_region4[9]+n_region4[10]*D)),0.50))/2.00;
    t_Output = temperature;

    qDebug() << temperature;
}

void WaterSteam::reg_4_saturation_Pressure()
{
    double t_Ref = t_Reference_region4;
    double p_Ref = p_Reference_region4;
    double temperature = t_Input;
    double tetha = (temperature/t_Ref)+(n_region4[9]/((temperature/t_Ref)-n_region4[10]));

    double A = pow(tetha,2)+n_region4[1]*tetha+n_region4[2];
    double B = n_region4[3]*pow(tetha,2)+n_region4[4]*tetha+n_region4[5];
    double C = n_region4[6]*pow(tetha,2)+n_region4[7]*tetha+n_region4[8];

    double pressure = p_Ref*pow(((2.00*C)/(-B+sqrt(pow(B,2)-4.00*A*C))),4);
    p_Output = pressure;
}



void WaterSteam::init()
{
    ui->TemperatureIncomboBox->addItems(TemperatureUnit);
    ui->PressureIncomboBox->addItems(PressureUnit);
    ui->TemperatureOutcomboBox->addItems(TemperatureUnit);
    ui->PressureOutcomboBox->addItems(PressureUnit);
    ui->TemperatureSpinBox->setValue(0);
    ui->PressureSpinBox->setValue(0);
    ui->TemperaturelineEdit->setText("");
    ui->PressurelineEdit->setText("");
    ui->TemperatureIncomboBox->setCurrentIndex(0);
    ui->PressureIncomboBox->setCurrentIndex(0);
    ui->TemperatureOutcomboBox->setCurrentIndex(0);
    ui->PressureOutcomboBox->setCurrentIndex(0);
    ui->Notiflabel->setText("");
}

void WaterSteam::read_InputValue()
{
    double temperature = ui->TemperatureSpinBox->value();
    double pressure = ui->PressureSpinBox->value();

    int t = ui->TemperatureIncomboBox->currentIndex();
    int p = ui->PressureIncomboBox->currentIndex();

    if (t<0) {
        t = 0;
    } else {
        t = t;
    }

    if (p<0) {
        p = 0;
    } else {
        p = p;
    }

    t_Input = t_Conversion_Any_to_K(t, temperature);
    p_Input = pressure*PressureUnitCoefficient.at(p);
}

void WaterSteam::send_OutputValue()
{
    double temperature = t_Output;
    double pressure = p_Output;

    int t = ui->TemperatureOutcomboBox->currentIndex();
    int p = ui->PressureOutcomboBox->currentIndex();

    if (t<0) {
        t = 0;
    } else {
        t = t;
    }

    if (p<0) {
        p = 0;
    } else {
        p = p;
    }

    t_Output = t_Conversion_K_to_Any(t, temperature);
    p_Output = pressure/PressureUnitCoefficient.at(p);
}

void WaterSteam::print_TemperatureOutputValue()
{
    QString temp = QString::number(t_Output, 'g', 6);
    ui->TemperaturelineEdit->setText(temp);
}

void WaterSteam::print_PressureOutputValue()
{
    QString press = QString::number(p_Output, 'g', 6);
    ui->PressurelineEdit->setText(press);
}

void WaterSteam::update_calculation()
{
    read_InputValue();
    int i = ui->TemperatureOutcomboBox->currentIndex();

    if (i<0) {
        i = 0;
    } else {
        i = i;
    }

    double temperature = t_Conversion_K_to_Any(i, t_Input);
    QString temp = QString::number(temperature, 'g', 6);
    ui->TemperaturelineEdit->setText(temp);

    int j = ui->PressureOutcomboBox->currentIndex();


    if (j<0) {
        j = 0;
    } else {
        j = j;
    }

    double pressure = p_Input/PressureUnitCoefficient.at(j);
    QString press = QString::number(pressure, 'g', 6);
    ui->PressurelineEdit->setText(press);

    if (p_Input==0 && t_Input!=0) {
        reg_4_saturation_Pressure();
        send_OutputValue();
        print_PressureOutputValue();
    } else if (t_Input==273.15 && p_Input!=0) {
        reg_4_saturation_Temperature();
        send_OutputValue();
        print_TemperatureOutputValue();
        qDebug() << "ABS";
    }

    qDebug() << t_Input << p_Input;
}

void WaterSteam::on_TemperatureSpinBox_valueChanged()
{
    update_calculation();
}

void WaterSteam::on_PressureSpinBox_valueChanged()
{
    update_calculation();
}

void WaterSteam::on_TemperatureIncomboBox_currentIndexChanged()
{
    update_calculation();
}

void WaterSteam::on_PressureIncomboBox_currentIndexChanged()
{
    update_calculation();
}

void WaterSteam::on_TemperatureOutcomboBox_currentIndexChanged()
{
    update_calculation();
}

void WaterSteam::on_PressureOutcomboBox_currentIndexChanged()
{
    update_calculation();
}

void WaterSteam::on_resetButton_clicked()
{
    init();
}
