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
