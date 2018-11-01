#ifndef WATERSTEAM_H
#define WATERSTEAM_H

#include <QMainWindow>
#include "qmath.h"

namespace Ui {
class WaterSteam;
}

class WaterSteam : public QMainWindow
{
    Q_OBJECT

public:
    explicit WaterSteam(QWidget *parent = 0);
    ~WaterSteam();

    //REFERENCE
    double t_Reference_region4 = 1.00;
    double p_Reference_region4 = 1.00;

    //INPUT
    double t_Input;
    double p_Input;

    //OUTPUT
    double t_Output;
    double p_Output;

    //FUNCTION
    void reg_4_saturation_Temperature();
    void reg_4_saturation_Pressure();

    //CONVERSION
    //                                          #1          #2          #3              #4              #5              #6              #7              #8              #9              #10
    QList<QString> PressureUnit             =   {"MPa",     "kPa",      "Pa",           "mm H2O",       "in H2O",       "ft H2O",       "mm Hg",        "in Hg",        "psi",          "ksi"};
    QList<double> PressureUnitCoefficient   =   {1.00,      0.001,      0.000001,       0.000009,       0.00025,        0.00299,        0.00013,        0.00339,        0.00689,        6.89476};
    QList<double> PressureMin               =   {0.0006,    0.60,       600.00,         61.18466,       2.40885,        0.20074,        4.50038,        0.17718,        0.08702,        0.00008};
    QList<double> PressureMax               =   {100.00,    100000.00,  100000000.00,   10197442.889,   401474.21331,   33456.22921,    750063.75542,   29530.05865,    14503.77378,    14.50377};

    //                                          #1          #2              #3          #4
    QList<QString> TemperatureUnit          =   {"Celcius", "Fahrenheit",   "Kelvin",   "Rankine"};
    QList<double> TemperatureMin            =   {0.01,      32.018,         273.16,     491.688};
    QList<double> TemperatureMax            =   {800.00,    1472.00,        1073.15,    1931.67};



    //COEFFICIENT
    double n_region4[11] =   {
                        0.00,
                        1167.0521452767,
                        -724213.16703206,
                        -17.073846940092,
                        12020.82470247,
                        -3232555.0322333,
                        14.91510861353,
                        -4823.2657361591,
                        405113.40542057,
                        -0.23855557567849,
                        650.17534844798,
                        };

    //FUNCTION
    void init();
    void read_InputValue();
    void send_OutputValue();
    void print_TemperatureOutputValue();
    void print_PressureOutputValue();
    void update_calculation();


private slots:
    void on_TemperatureSpinBox_valueChanged();

    void on_PressureSpinBox_valueChanged();

    void on_TemperatureIncomboBox_currentIndexChanged();

    void on_PressureIncomboBox_currentIndexChanged();

    void on_TemperatureOutcomboBox_currentIndexChanged();

    void on_PressureOutcomboBox_currentIndexChanged();

    void on_resetButton_clicked();

private:
    Ui::WaterSteam *ui;
};

#endif // WATERSTEAM_H
