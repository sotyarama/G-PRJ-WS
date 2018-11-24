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



    //QComboBox* aCombo = new QComboBox;
    //aCombo->addItems(TemperatureUnit);

    //reg_2_metastable_thermodynamic_Properties();
    //reg_1_thermodynamic_Properties();
    //print_region1_table();
    //ui->tableWidget->setColumnCount(2);
    //ui->tableWidget->setRowCount(2);
    //ui->tableWidget->setCellWidget(1,1, aCombo);
}

WaterSteam::~WaterSteam()
{
    delete ui;
}

//[REGION - 1]
void WaterSteam::reg_1_thermodynamic_Properties()
{
    double t_Ref = 1386.00;
    double p_Ref = 16.53;
    double tau = t_Ref/t_Input1;
    double pi = p_Input1/p_Ref;

    double gamma = 0.00;
    double gamma_p = 0.00;
    double gamma_pp = 0.00;
    double gamma_t = 0.00;
    double gamma_tt = 0.00;
    double gamma_p_t = 0.00;

    for (int i = 1; i < 35; i++) {
        double g = n_region1[i]*pow((7.1-pi),i_region1[i])*pow((tau-1.222),j_region1[i]);
        double g_p = -n_region1[i]*i_region1[i]*pow((7.1-pi),(i_region1[i]-1))*pow((tau-1.222),j_region1[i]);
        double g_pp = n_region1[i]*i_region1[i]*(i_region1[i]-1)*pow((7.1-pi),(i_region1[i]-2))*pow((tau-1.222),j_region1[i]);
        double g_t = n_region1[i]*j_region1[i]*pow((7.1-pi),i_region1[i])*pow((tau-1.222),(j_region1[i]-1));
        double g_tt = n_region1[i]*j_region1[i]*(j_region1[i]-1)*pow((7.1-pi),i_region1[i])*pow((tau-1.222),(j_region1[i]-2));
        double g_p_t = -n_region1[i]*i_region1[i]*j_region1[i]*pow((7.1-pi),(i_region1[i]-1))*pow((tau-1.222),(j_region1[i]-1));

        gamma = gamma+g;
        gamma_p = gamma_p+g_p;
        gamma_pp = gamma_pp+g_pp;
        gamma_t = gamma_t+g_t;
        gamma_tt = gamma_tt+g_tt;
        gamma_p_t = gamma_p_t+g_p_t;
    }

    v_1 = (pi*gamma_p*((r_Reference*t_Input1)/p_Input1))/1000;
    rho_1 = (p_Input1/(r_Reference*t_Input1*(pi*gamma_p)))*1000;
    u_1 = r_Reference*t_Input1*(tau*gamma_t-pi*gamma_p);
    s_1 = r_Reference*(tau*gamma_t-gamma);
    h_1 = r_Reference*t_Input1*tau*gamma_t;
    cv_1 = r_Reference*((-pow(tau,2))*gamma_tt+(pow((gamma_p-tau*gamma_p_t),2)/gamma_pp));
    cp_1 = r_Reference*-pow(tau,2)*gamma_tt;
}

//[REGION - 2]
void WaterSteam::reg_2_thermodynamic_Properties()
{
    double t_Ref = 540.00;
    double p_Ref = 1.00;
    double tau = t_Ref/t_Input2;
    double pi = p_Input2/p_Ref;

    double i_gamma = 0.00;
    double i_gamma_p = 0.00;
    double i_gamma_pp = 0.00;
    double i_gamma_t = 0.00;
    double i_gamma_tt = 0.00;
    double i_gamma_p_t = 0.00;

    double r_gamma = 0.00;
    double r_gamma_p = 0.00;
    double r_gamma_pp = 0.00;
    double r_gamma_t = 0.00;
    double r_gamma_tt = 0.00;
    double r_gamma_p_t = 0.00;

    for (int i = 1; i < 10; i++) {
        double g = n0_region2[i]*pow(tau,j0_region2[i]);
        double g_t = n0_region2[i]*j0_region2[i]*pow(tau,(j0_region2[i]-1));
        double g_tt = n0_region2[i]*j0_region2[i]*(j0_region2[i]-1)*pow(tau,(j0_region2[i]-2));

        i_gamma = i_gamma+g;
        i_gamma_t = i_gamma_t+g_t;
        i_gamma_tt = i_gamma_tt+g_tt;
    }

    i_gamma = i_gamma+log(pi);
    i_gamma_p = 1/pi;
    i_gamma_pp = -1/pow(pi,2);

    for (int i = 1; i < 44; i++) {
        double g = n_region2[i]*pow(pi,i_region2[i])*pow((tau-0.5),j_region2[i]);
        double g_p = n_region2[i]*i_region2[i]*pow(pi,(i_region2[i]-1))*pow((tau-0.5),j_region2[i]);
        double g_pp = n_region2[i]*i_region2[i]*(i_region2[i]-1)*pow(pi,(i_region2[i]-2))*pow((tau-0.5),j_region2[i]);
        double g_t = n_region2[i]*pow(pi,i_region2[i])*j_region2[i]*pow((tau-0.5),(j_region2[i]-1));
        double g_tt = n_region2[i]*pow(pi,i_region2[i])*j_region2[i]*(j_region2[i]-1)*pow((tau-0.5),(j_region2[i]-2));
        double g_p_t = n_region2[i]*i_region2[i]*pow(pi,(i_region2[i]-1))*j_region2[i]*pow((tau-0.5),(j_region2[i]-1));

        r_gamma = r_gamma+g;
        r_gamma_p = r_gamma_p+g_p;
        r_gamma_pp = r_gamma_pp+g_pp;
        r_gamma_t = r_gamma_t+g_t;
        r_gamma_tt = r_gamma_tt+g_tt;
        r_gamma_p_t = r_gamma_p_t+g_p_t;
    }

    v_2 = (pi*(i_gamma_p+r_gamma_p)*((r_Reference*t_Input2)/p_Input2))/1000;
    rho_2 = (p_Input2/(r_Reference*t_Input2*pi*(i_gamma_p+r_gamma_p)))*1000;
    u_2 = r_Reference*t_Input2*(tau*(i_gamma_t+r_gamma_t)-pi*(i_gamma_p+r_gamma_p));
    s_2 = r_Reference*(tau*(i_gamma_t+r_gamma_t)-(i_gamma+r_gamma));
    h_2 = r_Reference*t_Input2*tau*(i_gamma_t+r_gamma_t);
    cv_2 = r_Reference*(-pow(tau,2)*(i_gamma_tt+r_gamma_tt)-(pow((1+pi*r_gamma_p-tau*pi*r_gamma_p_t),2)/(1-pow(pi,2)*r_gamma_pp)));
    cp_2 = r_Reference*-pow(tau,2)*(i_gamma_tt+r_gamma_tt);
}

//[REGION - 2 Metastable]
void WaterSteam::reg_2_metastable_thermodynamic_Properties()
{
    double t_Ref = 540;
    double p_Ref = 1.00;
    double tau = t_Ref/t_Input2meta;
    double pi = p_Input2meta/p_Ref;

    double i_gamma = 0.00;
    double i_gamma_p = 0.00;
    double i_gamma_pp = 0.00;
    double i_gamma_t = 0.00;
    double i_gamma_tt = 0.00;
    double i_gamma_p_t = 0.00;

    double r_gamma = 0.00;
    double r_gamma_p = 0.00;
    double r_gamma_pp = 0.00;
    double r_gamma_t = 0.00;
    double r_gamma_tt = 0.00;
    double r_gamma_p_t = 0.00;

    for (int i = 1; i < 10; i++) {
        double g = n0_region2meta[i]*pow(tau,j0_region2meta[i]);
        double g_t = n0_region2meta[i]*j0_region2meta[i]*pow(tau,(j0_region2meta[i]-1));
        double g_tt = n0_region2meta[i]*j0_region2meta[i]*(j0_region2meta[i]-1)*pow(tau,(j0_region2meta[i]-2));

        i_gamma = i_gamma+g;
        i_gamma_t = i_gamma_t+g_t;
        i_gamma_tt = i_gamma_tt+g_tt;
    }

    i_gamma = log(pi)+i_gamma;
    i_gamma_p = 1/pi;
    i_gamma_pp = -1/pow(pi,2);

    for (int i = 1; i < 14; i++) {
        double g = n_region2meta[i]*pow(pi,i_region2meta[i])*pow((tau-0.5),j_region2meta[i]);
        double g_p = n_region2meta[i]*i_region2meta[i]*pow(pi,(i_region2meta[i]-1))*pow((tau-0.5),j_region2meta[i]);
        double g_pp = n_region2meta[i]*i_region2meta[i]*(i_region2meta[i]-1)*pow(pi,(i_region2meta[i]-2))*pow((tau-0.5),j_region2meta[i]);
        double g_t = n_region2meta[i]*pow(pi,i_region2meta[i])*j_region2meta[i]*pow((tau-0.5),(j_region2meta[i]-1));
        double g_tt = n_region2meta[i]*pow(pi,i_region2meta[i])*j_region2meta[i]*(j_region2meta[i]-1)*pow((tau-0.5),(j_region2meta[i]-2));
        double g_p_t = n_region2meta[i]*i_region2meta[i]*pow(pi,(i_region2meta[i]-1))*j_region2meta[i]*pow((tau-0.5),(j_region2meta[i]-1));

        r_gamma = r_gamma+g;
        r_gamma_p = r_gamma_p+g_p;
        r_gamma_pp = r_gamma_pp+g_pp;
        r_gamma_t = r_gamma_t+g_t;
        r_gamma_tt = r_gamma_tt+g_tt;
        r_gamma_p_t = r_gamma_p_t+g_p_t;
    }

    v_2meta = (pi*(i_gamma_p+r_gamma_p)*((r_Reference*t_Input2meta)/p_Input2meta))/1000;
    rho_2meta = (p_Input2meta/(r_Reference*t_Input2meta*(pi*(i_gamma_p+r_gamma_p))))*1000;
    u_2meta = r_Reference*t_Input2meta*(tau*(i_gamma_t+r_gamma_t)-pi*(i_gamma_p+r_gamma_p));
    s_2meta = r_Reference*(tau*(i_gamma_t+r_gamma_t)-(i_gamma+r_gamma));
    h_2meta = r_Reference*t_Input2meta*tau*(i_gamma_t+r_gamma_t);
    cv_2meta = r_Reference*(-pow(tau,2)*(i_gamma_tt+r_gamma_tt)-(pow((1+pi*r_gamma_p-tau*pi*r_gamma_p_t),2)/(1-pow(pi,2)*r_gamma_pp)));
    cp_2meta = r_Reference*(-pow(tau,2)*(i_gamma_tt+r_gamma_tt));
}

//[REGION - 4]
void WaterSteam::reg_4_saturation()
{
    double t_Ref = 1.00;
    double p_Ref = 1.00;
    double temperature = t_Input;
    double pressure = p_Input;
    double beta = pow((pressure/p_Ref),0.25);
    double tetha = (temperature/t_Ref)+(n_region4[9]/((temperature/t_Ref)-n_region4[10]));

    double A = pow(tetha,2)+n_region4[1]*tetha+n_region4[2];
    double B = n_region4[3]*pow(tetha,2)+n_region4[4]*tetha+n_region4[5];
    double C = n_region4[6]*pow(tetha,2)+n_region4[7]*tetha+n_region4[8];
    double E = pow(beta,2)+n_region4[3]*beta+n_region4[6];
    double F = n_region4[1]*pow(beta,2)+n_region4[4]*beta+n_region4[7];
    double G = n_region4[2]*pow(beta,2)+n_region4[5]*beta+n_region4[8];
    double D = (2.00*G)/(-F-pow((pow(F,2)-4.00*E*G),0.50));

    double t = t_Ref*(n_region4[10]+D-pow((pow((n_region4[10]+D),2)-4.00*(n_region4[9]+n_region4[10]*D)),0.50))/2.00;
    double p = p_Ref*pow(((2.00*C)/(-B+sqrt(pow(B,2)-4.00*A*C))),4);
    t_Output = t;
    p_Output = p;
}

//[REGION - 5]
void WaterSteam::reg_5_thermodynamic_Properties()
{
    double t_Ref = 1000.00;
    double p_Ref = 1.00;
    double tau = t_Ref/t_Input5;
    double pi = p_Input5/p_Ref;

    double i_gamma = 0.00;
    double i_gamma_p = 0.00;
    double i_gamma_pp = 0.00;
    double i_gamma_t = 0.00;
    double i_gamma_tt = 0.00;
    double i_gamma_p_t = 0.00;

    double r_gamma = 0.00;
    double r_gamma_p = 0.00;
    double r_gamma_pp = 0.00;
    double r_gamma_t = 0.00;
    double r_gamma_tt = 0.00;
    double r_gamma_p_t = 0.00;

    for (int i = 1; i < 7; i++) {
        double i_g = n0_region5[i]*pow(tau,j0_region5[i]);
        double i_g_t = n0_region5[i]*j0_region5[i]*pow(tau,(j0_region5[i]-1));
        double i_g_tt = n0_region5[i]*j0_region5[i]*(j0_region5[i]-1)*pow(tau,(j0_region5[i]-2));

        i_gamma = i_gamma+i_g;
        i_gamma_t = i_gamma_t+i_g_t;
        i_gamma_tt = i_gamma_tt+i_g_tt;

        double r_g = n_region5[i]*pow(pi,i_region5[i])*pow(tau,j_region5[i]);
        double r_g_p = n_region5[i]*i_region5[i]*pow(pi,(i_region5[i]-1))*pow(tau,j_region5[i]);
        double r_g_pp = n_region5[i]*i_region5[i]*(i_region5[i]-1)*pow(pi,(i_region5[i]-2))*pow(tau,j_region5[i]);
        double r_g_t = n_region5[i]*pow(pi,i_region5[i])*j_region5[i]*pow(tau,(j_region5[i]-1));
        double r_g_tt = n_region5[i]*pow(pi,i_region5[i])*j_region5[i]*(j_region5[i]-1)*pow(tau,(j_region5[i]-2));
        double r_g_p_t = n_region5[i]*i_region5[i]*pow(pi,(i_region5[i]-1))*j_region5[i]*pow(tau,(j_region5[i]-1));

        r_gamma = r_gamma+r_g;
        r_gamma_p = r_gamma_p+r_g_p;
        r_gamma_pp = r_gamma_pp+r_g_pp;
        r_gamma_t = r_gamma_t+r_g_t;
        r_gamma_tt = r_gamma_tt+r_g_tt;
        r_gamma_p_t = r_gamma_p_t+r_g_p_t;
    }

    i_gamma = log(pi)+i_gamma;
    i_gamma_p = 1/pi;
    i_gamma_pp = -1/pow(pi,2);
}

//[GENERAL - FUNCTION]

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
    ui->Notiflabel->setText("Please enter either temperature or pressure or both");
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

void WaterSteam::print_TemperatureInputValue()
{
    int i = ui->TemperatureOutcomboBox->currentIndex();

    if (i<0) {
        i = 0;
    } else {
        i = i;
    }

    double temperature = t_Conversion_K_to_Any(i, t_Input);
    QString temp = QString::number(temperature, 'g', 6);
    ui->TemperaturelineEdit->setText(temp);
}

void WaterSteam::print_PressureInputValue()
{
    int i = ui->PressureOutcomboBox->currentIndex();

    if (i<0) {
        i = 0;
    } else {
        i = i;
    }

    double pressure = p_Input/PressureUnitCoefficient.at(i);
    QString press = QString::number(pressure, 'g', 6);
    ui->PressurelineEdit->setText(press);
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

void WaterSteam::print_region1_table()
{
    QList<QString> tableH_Header = {"Properties", "Liquid", "Unit"};
    QString styleSheet = "::section {"
                         "background-color: blue;}";

    double value[7] = {
                    h_1, s_1, u_1, rho_1, v_1, cv_1, cp_1
                    };
    QStringList label = {"Enthalpy", "Entropy", "Internal Energy", "Density", "Volume", "Isochoric Heat Capacity", "Isobaric Heat Capacity"};

    ui->tableWidget->setRowCount(7);
    ui->tableWidget->setColumnCount(3);
    ui->tableWidget->setHorizontalHeaderLabels(tableH_Header);
    ui->tableWidget->horizontalHeader()->setStyleSheet(styleSheet);

    for (int i = 0; i < 7; i++) {
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(label[i]));
        ui->tableWidget->setItem(i, 1, new QTableWidgetItem(QString::number(value[i], 'e', 9)));
    }

    QComboBox* h = new QComboBox;
    QComboBox* s = new QComboBox;
    QComboBox* u = new QComboBox;
    QComboBox* rho = new QComboBox;
    QComboBox* v = new QComboBox;
    QComboBox* cv = new QComboBox;
    QComboBox * cp = new QComboBox;

    h->addItems(EnthalpyUnit);
    s->addItems(EntropyUnit);
    u->addItems(EnthalpyUnit);
    rho->addItems(DensityUnit);
    v->addItems(VolumeUnit);
    cv->addItems(EntropyUnit);
    cp->addItems(EntropyUnit);

    ui->tableWidget->setMinimumWidth(420);

    ui->tableWidget->setColumnWidth(0, 160);
    ui->tableWidget->setColumnWidth(1, 130);
    ui->tableWidget->setColumnWidth(2, 130);

    ui->tableWidget->setCellWidget(0, 2, h);
    ui->tableWidget->setCellWidget(1, 2, s);
    ui->tableWidget->setCellWidget(2, 2, u);
    ui->tableWidget->setCellWidget(3, 2, rho);
    ui->tableWidget->setCellWidget(4, 2, v);
    ui->tableWidget->setCellWidget(5, 2, cv);
    ui->tableWidget->setCellWidget(6, 2, cp);

    connect(h, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double H = h_1*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(0, 1, new QTableWidgetItem(QString::number(H, 'g', 12)));
    });
    connect(s, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double S = s_1*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(1, 1, new QTableWidgetItem(QString::number(S, 'g', 12)));
    });
    connect(u, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double U = u_1*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(2, 1, new QTableWidgetItem(QString::number(U, 'g', 12)));
    });
    connect(rho, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double RHO = rho_1*DensityUnitCoefficient.at(index);
        ui->tableWidget->setItem(3, 1, new QTableWidgetItem(QString::number(RHO, 'g', 12)));
    });
    connect(v, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double V = v_1*VolumeUnitCoefficient.at(index);
        ui->tableWidget->setItem(4, 1, new QTableWidgetItem(QString::number(V, 'g', 12)));
    });
    connect(cv, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CV = cv_1*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(5, 1, new QTableWidgetItem(QString::number(CV, 'g', 12)));
    });
    connect(cp, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CP = cp_1*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(6, 1, new QTableWidgetItem(QString::number(CP, 'g', 12)));
    });

}

void WaterSteam::print_region2_table(){
    QList<QString> tableH_Header = {"Properties", "Vapor", "Unit"};
    QString styleSheet = "::section {"
                         "background-color: blue;}";

    double value[7] = {
                    h_2, s_2, u_2, rho_2, v_2, cv_2, cp_2
                    };

    QStringList label = {"Enthalpy", "Entropy", "Internal Energy", "Density", "Volume", "Isochoric Heat Capacity", "Isobaric Heat Capacity"};

    ui->tableWidget->setRowCount(7);
    ui->tableWidget->setColumnCount(3);
    ui->tableWidget->setHorizontalHeaderLabels(tableH_Header);
    ui->tableWidget->horizontalHeader()->setStyleSheet(styleSheet);

    for (int i = 0; i < 7; i++) {
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(label[i]));
        ui->tableWidget->setItem(i, 1, new QTableWidgetItem(QString::number(value[i], 'e', 9)));
    }

    QComboBox* h = new QComboBox;
    QComboBox* s = new QComboBox;
    QComboBox* u = new QComboBox;
    QComboBox* rho = new QComboBox;
    QComboBox* v = new QComboBox;
    QComboBox* cv = new QComboBox;
    QComboBox * cp = new QComboBox;

    h->addItems(EnthalpyUnit);
    s->addItems(EntropyUnit);
    u->addItems(EnthalpyUnit);
    rho->addItems(DensityUnit);
    v->addItems(VolumeUnit);
    cv->addItems(EntropyUnit);
    cp->addItems(EntropyUnit);

    ui->tableWidget->setMinimumWidth(420);

    ui->tableWidget->setColumnWidth(0, 160);
    ui->tableWidget->setColumnWidth(1, 130);
    ui->tableWidget->setColumnWidth(2, 130);

    ui->tableWidget->setCellWidget(0, 2, h);
    ui->tableWidget->setCellWidget(1, 2, s);
    ui->tableWidget->setCellWidget(2, 2, u);
    ui->tableWidget->setCellWidget(3, 2, rho);
    ui->tableWidget->setCellWidget(4, 2, v);
    ui->tableWidget->setCellWidget(5, 2, cv);
    ui->tableWidget->setCellWidget(6, 2, cp);

    connect(h, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double H = h_2*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(0, 1, new QTableWidgetItem(QString::number(H, 'g', 12)));
    });
    connect(s, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double S = s_2*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(1, 1, new QTableWidgetItem(QString::number(S, 'g', 12)));
    });
    connect(u, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double U = u_2*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(2, 1, new QTableWidgetItem(QString::number(U, 'g', 12)));
    });
    connect(rho, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double RHO = rho_2*DensityUnitCoefficient.at(index);
        ui->tableWidget->setItem(3, 1, new QTableWidgetItem(QString::number(RHO, 'g', 12)));
    });
    connect(v, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double V = v_2*VolumeUnitCoefficient.at(index);
        ui->tableWidget->setItem(4, 1, new QTableWidgetItem(QString::number(V, 'g', 12)));
    });
    connect(cv, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CV = cv_2*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(5, 1, new QTableWidgetItem(QString::number(CV, 'g', 12)));
    });
    connect(cp, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CP = cp_2*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(6, 1, new QTableWidgetItem(QString::number(CP, 'g', 12)));
    });
}

void WaterSteam::print_region2meta_table()
{
    QList<QString> tableH_Header = {"Properties", "Vapor", "Unit"};
    QString styleSheet = "::section {"
                         "background-color: blue;}";

    double value[7] = {
                    h_2meta, s_2meta, u_2meta, rho_2meta, v_2meta, cv_2meta, cp_2meta
                    };
    QStringList label = {"Enthalpy", "Entropy", "Internal Energy", "Density", "Volume", "Isochoric Heat Capacity", "Isobaric Heat Capacity"};

    ui->tableWidget->setRowCount(7);
    ui->tableWidget->setColumnCount(3);
    ui->tableWidget->setHorizontalHeaderLabels(tableH_Header);
    ui->tableWidget->horizontalHeader()->setStyleSheet(styleSheet);

    for (int i = 0; i < 7; i++) {
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(label[i]));
        ui->tableWidget->setItem(i, 1, new QTableWidgetItem(QString::number(value[i], 'e', 9)));
    }

    QComboBox* h = new QComboBox;
    QComboBox* s = new QComboBox;
    QComboBox* u = new QComboBox;
    QComboBox* rho = new QComboBox;
    QComboBox* v = new QComboBox;
    QComboBox* cv = new QComboBox;
    QComboBox * cp = new QComboBox;

    h->addItems(EnthalpyUnit);
    s->addItems(EntropyUnit);
    u->addItems(EnthalpyUnit);
    rho->addItems(DensityUnit);
    v->addItems(VolumeUnit);
    cv->addItems(EntropyUnit);
    cp->addItems(EntropyUnit);

    ui->tableWidget->setMinimumWidth(420);

    ui->tableWidget->setColumnWidth(0, 160);
    ui->tableWidget->setColumnWidth(1, 130);
    ui->tableWidget->setColumnWidth(2, 130);

    ui->tableWidget->setCellWidget(0, 2, h);
    ui->tableWidget->setCellWidget(1, 2, s);
    ui->tableWidget->setCellWidget(2, 2, u);
    ui->tableWidget->setCellWidget(3, 2, rho);
    ui->tableWidget->setCellWidget(4, 2, v);
    ui->tableWidget->setCellWidget(5, 2, cv);
    ui->tableWidget->setCellWidget(6, 2, cp);

    connect(h, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double H = h_2meta*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(0, 1, new QTableWidgetItem(QString::number(H, 'e', 9)));
    });
    connect(s, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double S = s_2meta*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(1, 1, new QTableWidgetItem(QString::number(S, 'e', 9)));
    });
    connect(u, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double U = u_2meta*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(2, 1, new QTableWidgetItem(QString::number(U, 'e', 9)));
    });
    connect(rho, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double RHO = rho_2meta*DensityUnitCoefficient.at(index);
        ui->tableWidget->setItem(3, 1, new QTableWidgetItem(QString::number(RHO, 'e', 9)));
    });
    connect(v, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double V = v_2meta*VolumeUnitCoefficient.at(index);
        ui->tableWidget->setItem(4, 1, new QTableWidgetItem(QString::number(V, 'e', 9)));
    });
    connect(cv, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CV = cv_2meta*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(5, 1, new QTableWidgetItem(QString::number(CV, 'e', 9)));
    });
    connect(cp, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CP = cp_2meta*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(6, 1, new QTableWidgetItem(QString::number(CP, 'e', 9)));
    });
}

void WaterSteam::print_region4_table()
{
    QList<QString> tableH_Header = {"Properties", "Liquid", "Vapor", "Unit"};
    QString styleSheet = "::section {"
                         "background-color: blue;}";

    double value1[7] = {
                    h_1, s_1, u_1, rho_1, v_1, cv_1, cp_1
                    };
    double value2[7] = {
                    h_2meta, s_2meta, u_2meta, rho_2meta, v_2meta, cv_2meta, cp_2meta
                    };
    QStringList label = {"Enthalpy", "Entropy", "Internal Energy", "Density", "Volume", "Isochoric Heat Capacity", "Isobaric Heat Capacity"};

    ui->tableWidget->setRowCount(7);
    ui->tableWidget->setColumnCount(4);
    ui->tableWidget->setHorizontalHeaderLabels(tableH_Header);
    ui->tableWidget->horizontalHeader()->setStyleSheet(styleSheet);

    for (int i = 0; i < 7; i++) {
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(label[i]));
        ui->tableWidget->setItem(i, 1, new QTableWidgetItem(QString::number(value1[i], 'e', 9)));
        ui->tableWidget->setItem(i, 2, new QTableWidgetItem(QString::number(value2[i], 'e', 9)));
    }

    QComboBox* h = new QComboBox;
    QComboBox* s = new QComboBox;
    QComboBox* u = new QComboBox;
    QComboBox* rho = new QComboBox;
    QComboBox* v = new QComboBox;
    QComboBox* cv = new QComboBox;
    QComboBox * cp = new QComboBox;

    h->addItems(EnthalpyUnit);
    s->addItems(EntropyUnit);
    u->addItems(EnthalpyUnit);
    rho->addItems(DensityUnit);
    v->addItems(VolumeUnit);
    cv->addItems(EntropyUnit);
    cp->addItems(EntropyUnit);

    ui->tableWidget->setMinimumWidth(550);

    ui->tableWidget->setColumnWidth(0, 160);
    ui->tableWidget->setColumnWidth(1, 130);
    ui->tableWidget->setColumnWidth(2, 130);
    ui->tableWidget->setColumnWidth(3, 130);

    ui->tableWidget->setCellWidget(0, 3, h);
    ui->tableWidget->setCellWidget(1, 3, s);
    ui->tableWidget->setCellWidget(2, 3, u);
    ui->tableWidget->setCellWidget(3, 3, rho);
    ui->tableWidget->setCellWidget(4, 3, v);
    ui->tableWidget->setCellWidget(5, 3, cv);
    ui->tableWidget->setCellWidget(6, 3, cp);

    connect(h, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double H_1 = h_1*EnthalpyUnitCoefficient.at(index);
        double H_2meta = h_2meta*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(0, 1, new QTableWidgetItem(QString::number(H_1, 'e', 9)));
        ui->tableWidget->setItem(0, 2, new QTableWidgetItem(QString::number(H_2meta, 'e', 9)));
    });
    connect(s, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double S_1 = s_1*EntropyUnitCoefficient.at(index);
        double S_2meta = s_2meta*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(1, 1, new QTableWidgetItem(QString::number(S_1, 'e', 9)));
        ui->tableWidget->setItem(1, 2, new QTableWidgetItem(QString::number(S_2meta, 'e', 9)));
    });
    connect(u, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double U_1 = u_1*EnthalpyUnitCoefficient.at(index);
        double U_2meta = u_2meta*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(2, 1, new QTableWidgetItem(QString::number(U_1, 'e', 9)));
        ui->tableWidget->setItem(2, 2, new QTableWidgetItem(QString::number(U_2meta, 'e', 9)));
    });
    connect(rho, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double RHO_1 = rho_1*DensityUnitCoefficient.at(index);
        double RHO_2meta = rho_2meta*DensityUnitCoefficient.at(index);
        ui->tableWidget->setItem(3, 1, new QTableWidgetItem(QString::number(RHO_1, 'e', 9)));
        ui->tableWidget->setItem(3, 2, new QTableWidgetItem(QString::number(RHO_2meta, 'e', 9)));
    });
    connect(v, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double V_1 = v_1*VolumeUnitCoefficient.at(index);
        double V_2meta = v_2meta*VolumeUnitCoefficient.at(index);
        ui->tableWidget->setItem(4, 1, new QTableWidgetItem(QString::number(V_1, 'e', 9)));
        ui->tableWidget->setItem(4, 2, new QTableWidgetItem(QString::number(V_2meta, 'e', 9)));
    });
    connect(cv, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CV_1 = cv_1*EntropyUnitCoefficient.at(index);
        double CV_2meta = cv_2meta*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(5, 1, new QTableWidgetItem(QString::number(CV_1, 'e', 9)));
        ui->tableWidget->setItem(5, 2, new QTableWidgetItem(QString::number(CV_2meta, 'e', 9)));
    });
    connect(cp, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CP_1 = cp_1*EntropyUnitCoefficient.at(index);
        double CP_2meta = cp_2meta*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(6, 1, new QTableWidgetItem(QString::number(CP_1, 'e', 9)));
        ui->tableWidget->setItem(6, 2, new QTableWidgetItem(QString::number(CP_2meta, 'e', 9)));
    });
}

void WaterSteam::print_region5_table()
{
    QList<QString> tableH_Header = {"Properties", "Vapor", "Unit"};
    QString styleSheet = "::section {"
                         "background-color: blue;}";

    double value[7] = {
                    h_5, s_5, u_5, rho_5, v_5, cv_5, cp_5
                    };

    QStringList label = {"Enthalpy", "Entropy", "Internal Energy", "Density", "Volume", "Isochoric Heat Capacity", "Isobaric Heat Capacity"};

    ui->tableWidget->setRowCount(7);
    ui->tableWidget->setColumnCount(3);
    ui->tableWidget->setHorizontalHeaderLabels(tableH_Header);
    ui->tableWidget->horizontalHeader()->setStyleSheet(styleSheet);

    for (int i = 0; i < 7; i++) {
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(label[i]));
        ui->tableWidget->setItem(i, 1, new QTableWidgetItem(QString::number(value[i], 'e', 9)));
    }

    QComboBox* h = new QComboBox;
    QComboBox* s = new QComboBox;
    QComboBox* u = new QComboBox;
    QComboBox* rho = new QComboBox;
    QComboBox* v = new QComboBox;
    QComboBox* cv = new QComboBox;
    QComboBox * cp = new QComboBox;

    h->addItems(EnthalpyUnit);
    s->addItems(EntropyUnit);
    u->addItems(EnthalpyUnit);
    rho->addItems(DensityUnit);
    v->addItems(VolumeUnit);
    cv->addItems(EntropyUnit);
    cp->addItems(EntropyUnit);

    ui->tableWidget->setMinimumWidth(420);

    ui->tableWidget->setColumnWidth(0, 160);
    ui->tableWidget->setColumnWidth(1, 130);
    ui->tableWidget->setColumnWidth(2, 130);

    ui->tableWidget->setCellWidget(0, 2, h);
    ui->tableWidget->setCellWidget(1, 2, s);
    ui->tableWidget->setCellWidget(2, 2, u);
    ui->tableWidget->setCellWidget(3, 2, rho);
    ui->tableWidget->setCellWidget(4, 2, v);
    ui->tableWidget->setCellWidget(5, 2, cv);
    ui->tableWidget->setCellWidget(6, 2, cp);

    connect(h, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double H = h_5*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(0, 1, new QTableWidgetItem(QString::number(H, 'g', 12)));
    });
    connect(s, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double S = s_5*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(1, 1, new QTableWidgetItem(QString::number(S, 'g', 12)));
    });
    connect(u, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double U = u_5*EnthalpyUnitCoefficient.at(index);
        ui->tableWidget->setItem(2, 1, new QTableWidgetItem(QString::number(U, 'g', 12)));
    });
    connect(rho, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double RHO = rho_5*DensityUnitCoefficient.at(index);
        ui->tableWidget->setItem(3, 1, new QTableWidgetItem(QString::number(RHO, 'g', 12)));
    });
    connect(v, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double V = v_5*VolumeUnitCoefficient.at(index);
        ui->tableWidget->setItem(4, 1, new QTableWidgetItem(QString::number(V, 'g', 12)));
    });
    connect(cv, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CV = cv_5*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(5, 1, new QTableWidgetItem(QString::number(CV, 'g', 12)));
    });
    connect(cp, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        double CP = cp_5*EntropyUnitCoefficient.at(index);
        ui->tableWidget->setItem(6, 1, new QTableWidgetItem(QString::number(CP, 'g', 12)));
    });
}

void WaterSteam::update_calculation()
{
    read_InputValue();
    reg_4_saturation();
    send_OutputValue();

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
        print_PressureOutputValue();
    } else if (t_Input==273.15 && p_Input!=0) {
        print_TemperatureOutputValue();
        qDebug() << "ABS";
    }
}

void WaterSteam::on_TemperatureSpinBox_valueChanged(double arg)
{
    read_InputValue();

    print_TemperatureInputValue();

    double t = arg;
    double p = ui->PressureSpinBox->value();

    int t_index = ui->TemperatureIncomboBox->currentIndex();

    int tIndex = ui->TemperatureOutcomboBox->currentIndex();
    int pIndex = ui->PressureOutcomboBox->currentIndex();

    double t_InputMin = t_Conversion_K_to_Any(t_index, 273.15);
    double t_InputMax = t_Conversion_K_to_Any(t_index, 647.096);

    QString tMin = QString::number(t_InputMin, 'g', 6);
    QString tMax = QString::number(t_InputMax, 'g', 6);

    if (t==0.00 && p==0.00) {
        resetTable();
        ui->Notiflabel->setText("Please enter either temperature or pressure or both");
        ui->TemperaturelineEdit->setText("");
        ui->PressurelineEdit->setText("");
    } else if (t>=t_InputMin && t<=t_InputMax && p==0.00) {
        resetTable();
        reg_4_saturation();
        send_OutputValue();
        print_PressureOutputValue();

        t_Input1 = t_Input;
        t_Input2meta = t_Input;

        double pOutput = p_Output*PressureUnitCoefficient.at(pIndex);

        p_Input1 = pOutput;
        p_Input2meta = pOutput;

        reg_1_thermodynamic_Properties();
        reg_2_metastable_thermodynamic_Properties();

        print_region4_table();
        ui->Notiflabel->setText("REGION 4");
    } else if (t>=t_InputMin && t<=t_InputMax && p>0.00) {
        resetTable();
        reg_4_saturation();
        send_OutputValue();

        double pOutput = p_Output*PressureUnitCoefficient.at(pIndex);

        if (p_Input < pOutput) {
            resetTable();
            t_Input2 = t_Input;
            p_Input2 = p_Input;
            reg_2_thermodynamic_Properties();
            print_region2_table();
            ui->Notiflabel->setText("REGION 2");
        } else if (p_Input > pOutput) {
            resetTable();
            t_Input1 = t_Input;
            p_Input1 = p_Input;
            reg_1_thermodynamic_Properties();
            print_region1_table();
            ui->Notiflabel->setText("REGION 1");
        } else {
            qDebug() << "ANEEEEEEEEH";
        }
    } else {
        resetTable();
        ui->Notiflabel->setText("For saturation the temperature must be between " + tMin + " - " + tMax + "." );
        ui->TemperaturelineEdit->setText("");
        ui->PressurelineEdit->setText("");
    }
}

void WaterSteam::on_PressureSpinBox_valueChanged(double arg)
{
    read_InputValue();

    print_PressureInputValue();

    double p = arg;
    double t = ui->TemperatureSpinBox->value();

    int p_index = ui->PressureIncomboBox->currentIndex();

    int tIndex = ui->TemperatureOutcomboBox->currentIndex();
    int pIndex = ui->PressureOutcomboBox->currentIndex();

    double p_InputMin = 0.00*PressureUnitCoefficient.at(p_index);
    double p_InputMax = 100.00*PressureUnitCoefficient.at(p_index);

    QString pMin = QString::number(p_InputMin, 'g', 6);
    QString pMax = QString::number(p_InputMax, 'g', 6);

    if (t==0.00 && p==0.00) {
        resetTable();
        ui->Notiflabel->setText("Please enter either temperature or pressure or both");
        ui->TemperaturelineEdit->setText("");
        ui->PressurelineEdit->setText("");
    } else if (p>=p_InputMin && p<=p_InputMax && t==0.00) {
        resetTable();
        reg_4_saturation();
        send_OutputValue();
        print_TemperatureOutputValue();

        p_Input1 = p_Input;
        p_Input2meta = p_Input;

        double tOutput = t_Conversion_Any_to_K(tIndex, t_Output);

        t_Input1 = tOutput;
        t_Input2meta = tOutput;

        reg_1_thermodynamic_Properties();
        reg_2_metastable_thermodynamic_Properties();;

        print_region4_table();
        ui->Notiflabel->setText("REGION 4");
    } else if ((p>=p_InputMin || p<=p_InputMax) && t!=0.00) {
        resetTable();
        reg_4_saturation();
        send_OutputValue();

        double tOutput = t_Conversion_Any_to_K(tIndex, t_Output);

        if (t_Input < tOutput) {
            resetTable();
            t_Input1 = t_Input;
            p_Input1 = p_Input;
            reg_1_thermodynamic_Properties();
            print_region1_table();
            ui->Notiflabel->setText("REGION 1");
        } else if (t_Input > tOutput) {
            resetTable();
            t_Input2 = t_Input;
            p_Input2 = p_Input;
            reg_2_thermodynamic_Properties();
            print_region2_table();
            ui->Notiflabel->setText("REGION 2");
        } else {
            resetTable();
            reg_4_saturation();
            send_OutputValue();
            print_PressureOutputValue();

            t_Input1 = t_Input;
            t_Input2meta = t_Input;

            double pOutput = p_Output*PressureUnitCoefficient.at(pIndex);

            p_Input1 = pOutput;
            p_Input2meta = pOutput;

            reg_1_thermodynamic_Properties();
            reg_2_metastable_thermodynamic_Properties();

            print_region4_table();
            ui->Notiflabel->setText("REGION 4");

            qDebug() << "TAMBAH ANEH";
        }
    } else {
        resetTable();
        ui->Notiflabel->setText("For saturation the pressure must be between " + pMin + " - " + pMax + "." );
        ui->TemperaturelineEdit->setText("");
        ui->PressurelineEdit->setText("");
    }
}

void WaterSteam::on_TemperatureIncomboBox_currentIndexChanged()
{

}

void WaterSteam::on_PressureIncomboBox_currentIndexChanged()
{

}

void WaterSteam::on_TemperatureOutcomboBox_currentIndexChanged()
{

}

void WaterSteam::on_PressureOutcomboBox_currentIndexChanged()
{

}

void WaterSteam::on_resetButton_clicked()
{
    init();
}

void WaterSteam::resetTable()
{
    ui->tableWidget->clearContents();
    ui->tableWidget->setColumnCount(0);
    ui->tableWidget->setRowCount(0);
}
