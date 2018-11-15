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

    //LIMIT
    double t_low = 273.15;
    double t_mid = 1073.15;
    double t_high = 2273.15;

    double p_low = 0.00;
    double p_mid = 50.00;
    double p_high = 100.00;


    //REFERENCE
    double t_C = 647.096;
    double p_C = 22.064;
    double rho_C = 322.00;
    double r_Reference = 0.461526;

    //INPUT
    double t_Input;
    double t_Input1;
    double t_Input2;
    double t_Input2meta;
    double t_Input5;

    double p_Input;
    double p_Input1;
    double p_Input2;
    double p_Input2meta;
    double p_Input5;

    //OUTPUT
    double t_Output;
    double p_Output;

    double v_1;
    double rho_1;
    double u_1;
    double s_1;
    double h_1;
    double cv_1;
    double cp_1;

    double v_2;
    double rho_2;
    double u_2;
    double s_2;
    double h_2;
    double cv_2;
    double cp_2;

    double v_2meta;
    double rho_2meta;
    double u_2meta;
    double s_2meta;
    double h_2meta;
    double cv_2meta;
    double cp_2meta;

    double v_5;
    double rho_5;
    double u_5;
    double s_5;
    double h_5;
    double cv_5;
    double cp_5;

    //THERMODYNAMIC-FUNCTION
    void reg_1_thermodynamic_Properties();
    void reg_2_thermodynamic_Properties();
    void reg_2_metastable_thermodynamic_Properties();
    void reg_4_saturation();
    void reg_5_thermodynamic_Properties();

    //CONVERSION
    //                                              #1          #2          #3              #4              #5              #6              #7              #8              #9              #10
    QList<QString> PressureUnit                 =   {"MPa",     "kPa",      "Pa",           "mm H2O",       "in H2O",       "ft H2O",       "mm Hg",        "in Hg",        "psi",          "ksi"};
    QList<double> PressureUnitCoefficient       =   {1.00,      0.001,      0.000001,       0.000009,       0.00025,        0.00299,        0.00013,        0.00339,        0.00689,        6.89476};
    QList<double> PressureMin                   =   {0.0006,    0.60,       600.00,         61.18466,       2.40885,        0.20074,        4.50038,        0.17718,        0.08702,        0.00008};
    QList<double> PressureMax                   =   {100.00,    100000.00,  100000000.00,   10197442.889,   401474.21331,   33456.22921,    750063.75542,   29530.05865,    14503.77378,    14.50377};

    //                                              #1          #2              #3          #4
    QList<QString> TemperatureUnit              =   {"Celcius", "Fahrenheit",   "Kelvin",   "Rankine"};
    QList<double> TemperatureMin                =   {0.01,      32.018,         273.16,     491.688};
    QList<double> TemperatureMax                =   {800.00,    1472.00,        1073.15,    1931.67};

    //                                              #1          #2                  #3              #4          #5
    QList<QString> VolumeUnit                   =   {"m3/kg",   "liter/kg (cm3/g)", "US gal/lbm",   "in3/lbm",  "ft3/lbm"};
    QList<double> VolumeUnitCoefficient         =   {1.00,      1000.00,            119.8264,       27679.90,   16.01846};

    //
    QList<QString> DensityUnit                  =   {"kg/m3",   "kg/liter (g/cm3)", "lbm/US gal",   "lbm/in3",  "lbm/ft3"};
    QList<double> DensityUnitCoefficient        =   {1.00,      1000.00,            119.8264,       27679.90,   16.01846};

    //                                              #1          #2          #3                  #4              #5              #6              #7
    QList<QString> EnthalpyUnit                 =   {"kJ/kg",   "cal/g",    "psia/(lbm/ft3)",   "kW.h/lbm",     "hp.h/lbm",     "ft.lbf/lbm",   "Btu/lbm"};
    QList<double> EnthalpyUnitCoefficient       =   {1.00,      0.2388459,  2.323282,           0.0001259979,   0.0001689659,   334.5526,       0.4299226};

    //                                              #1              #2              #3                  #4                  #5              #6                  #7
    QList<QString> EntropyUnit                  =   {"kJ/(kg.K)",   "cal/(g.K)",    "bar.cm3/(g.K)",    "psia.ft3/(lbm.R)", "kW.h/(lbm.R)", "ft.lbf/(lbm.R)",   "Btu/(lbm.R)"};
    QList<double> EntropyUnitCoefficient        =   {1.00,          0.2388459,      10.00,              1.290712,           0.00006999882,  185.8625,           0.2388459};



    //COEFFICIENT
    double i_region1[35] =   {
                        0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
                        0.00,   0.00,   1.00,   1.00,   1.00,   1.00,   1.00,
                        1.00,   2.00,   2.00,   2.00,   2.00,   2.00,   3.00,
                        3.00,   3.00,   4.00,   4.00,   4.00,   5.00,   8.00,
                        8.00,   21.00,  23.00,  29.00,  30.00,  31.00,  32.00,
                        };

    double j_region1[35] =   {
                        0.00,   -2.00,  -1.00,  0.00,   1.00,   2.00,   3.00,
                        4.00,   5.00,   -9.00,  -7.00,  -1.00,  0.00,   1.00,
                        3.00,   -3.00,  0.00,   1.00,   3.00,   17.00,  -4.00,
                        0.00,   6.00,   -5.00,  -2.00,  10.00,  -8.00,  -11.00,
                        -6.00,  -29.00, -31.00, -38.00, -39.00, -40.00, -41.00,
                        };

    double n_region1[35] =   {
                        0.00,                   0.14632971213167,       -0.84548187169114,      -3.756360367204,        3.3855169168385,
                        -0.95791963387872,      0.15772038513228,       -0.016616417199501,     0.00081214629983568,    0.00028319080123804,
                        -0.00060706301565874,   -0.018990068218419,     -0.032529748770505,     -0.021841717175414,     -0.00005283835796993,
                        -0.00047184321073267,   -0.00030001780793026,   0.000047661393906987,   -4.4141845330846E-06,   -7.2694996297594E-16,
                        -0.000031679644845054,  -2.8270797985312E-06,   -8.5205128120103E-10,   -0.0000022425281908,    -6.5171222895601E-07,
                        -1.4341729937924E-13,   -4.0516996860117E-07,   -1.2734301741641E-09,   -1.7424871230634E-10,   -6.8762131295531E-19,
                        1.4478307828521E-20,    2.6335781662795E-23,    -1.1947622640071E-23,   1.8228094581404E-24,    -9.3537087292458E-26,
                        };

    double j0_region2[10] =  {
                        0.00,   0.00,   1.00,   -5.00,  -4.00,
                        -3.00,  -2.00,  -1.00,  2.00,   3.00,
                        };

    double n0_region2[10] =  {
                        0.00,               -9.6927686500217,	10.086655968018,	-0.005608791128302,	0.071452738081455,
                        -0.40710498223928,  1.4240819171444,	-4.383951131945,	-0.28408632460772,	0.021268463753307,
                        };

    double i_region2[44] =   {
                        0.00,	1.00,	1.00,	1.00,
                        1.00,	1.00,	2.00,	2.00,
                        2.00,	2.00,	2.00,	3.00,
                        3.00,	3.00,	3.00,	3.00,
                        4.00,	4.00,	4.00,	5.00,
                        6.00,	6.00,	6.00,	7.00,
                        7.00,	7.00,	8.00,	8.00,
                        9.00,	10.00,	10.00,	10.00,
                        16.00,	16.00,	18.00,	20.00,
                        20.00,	20.00,	21.00,	22.00,
                        23.00,	24.00,	24.00,	24.00,
                        };

    double j_region2[44] =   {
                        0.00,	0.00,	1.00,	2.00,
                        3.00,	6.00,	1.00,	2.00,
                        4.00,	7.00,	36.00,	0.00,
                        1.00,	3.00,	6.00,	35.00,
                        1.00,	2.00,	3.00,	7.00,
                        3.00,	16.00,	35.00,	0.00,
                        11.00,	25.00,	8.00,	36.00,
                        13.00,	4.00,	10.00,	14.00,
                        29.00,	50.00,	57.00,	20.00,
                        35.00,	48.00,	21.00,	53.00,
                        39.00,	26.00,	40.00,	58.00,
                        };

    double n_region2[44] =   {
                        0.00,                   -0.0017731742473213,	-0.017834862292358, 	-0.045996013696365,
                        -0.057581259083432, 	-0.05032527872793,      -0.000033032641670203,	-0.00018948987516315,
                        -0.0039392777243355,	-0.043797295650573,     -0.000026674547914087,	2.0481737692309E-08,
                        4.3870667284435E-07,	-0.00003227767723857,	-0.0015033924542148,	-0.040668253562649,
                        -7.8847309559367E-10,	1.2790717852285E-08,	4.8225372718507E-07,	2.2922076337661E-06,
                        -1.6714766451061E-11,	-0.0021171472321355,	-23.895741934104,       -5.905956432427E-18,
                        -1.2621808899101E-06,	-0.038946842435739,     1.1256211360459E-11,	-8.2311340897998,
                        1.9809712802088E-08,	1.0406965210174E-19,	-1.0234747095929E-13,	-1.0018179379511E-09,
                        -8.0882908646985E-11,	0.10693031879409,       -0.33662250574171,      8.9185845355421E-25,
                        3.0629316876232E-13,	-4.2002467698208E-06,	-5.9056029685639E-26,	3.7826947613457E-06,
                        -1.2768608934681E-15,	7.3087610595061E-29,	5.5414715350778E-17,	-9.436970724121E-07,
                        };

    double j0_region2meta[10] =   {
                            0.00,
                            0.00,
                            1.00,
                            -5.00,
                            -4.00,
                            -3.00,
                            -2.00,
                            -1.00,
                            2.00,
                            3.00,
                            };

    double n0_region2meta[10] =   {
                            0.00,
                            -9.6927686500217,
                            10.086655968018,
                            -0.005608791128302,
                            0.071452738081455,
                            -0.40710498223928,
                            1.4240819171444,
                            -4.383951131945,
                            -0.28408632460772,
                            0.021268463753307,
                            };

    double i_region2meta[14] =    {
                            0.00,
                            1.00,
                            1.00,
                            1.00,
                            1.00,
                            2.00,
                            2.00,
                            2.00,
                            3.00,
                            3.00,
                            4.00,
                            4.00,
                            5.00,
                            5.00,
                            };

    double j_region2meta[14] =    {
                            0.00,
                            0.00,
                            2.00,
                            5.00,
                            11.00,
                            1.00,
                            7.00,
                            16.00,
                            4.00,
                            16.00,
                            7.00,
                            10.00,
                            9.00,
                            10.00,
                            };

    double n_region2meta[14] =    {
                            0.00,
                            -0.0073362260186506,
                            -0.088223831943146,
                            -0.072334555213245,
                            -0.0040813178534455,
                            0.0020097803380207,
                            -0.053045921898642,
                            -0.007619040908697,
                            -0.0063498037657313,
                            -0.086043093028588,
                            0.007532158152277,
                            -0.0079238375446139,
                            -0.00022888160778447,
                            -0.002645650148281,
                            };

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

    double j0_region5[7] =   {
                        0.00,
                        0.00,
                        1.00,
                        -3.00,
                        -2.00,
                        -1.00,
                        2.00,
                        };

    double n0_region5[7] =   {
                        0.00,
                        -13.179983674201,
                        6.8540841634434,
                        -0.024805148933466,
                        0.36901534980333,
                        -3.1161318213925,
                        -0.32961626538917,
                        };

    double i_region5[7] =    {
                        0.00,
                        1.00,
                        1.00,
                        1.00,
                        2.00,
                        2.00,
                        3.00,
                        };

    double j_region5[7] =    {
                        0.00,
                        1.00,
                        2.00,
                        3.00,
                        3.00,
                        9.00,
                        7.00,
                        };

    double n_region5[7] =    {
                        0.00,
                        0.0015736404855259,
                        0.00090153761673944,
                        -0.0050270077677648,
                        2.2440037409485E-06,
                        -4.1163275453471E-06,
                        3.7919454822955E-08,
                        };

    //FUNCTION
    void init();
    void read_InputValue();
    void send_OutputValue();
    void print_TemperatureInputValue();
    void print_PressureInputValue();
    void print_TemperatureOutputValue();
    void print_PressureOutputValue();
    void print_region1_table();
    void print_region2_table();
    void print_region2meta_table();
    void print_region4_table();
    void print_region5_table();
    void update_calculation();
    void resetTable();


private slots:
    void on_TemperatureSpinBox_valueChanged(double arg);

    void on_PressureSpinBox_valueChanged(double arg);

    void on_TemperatureIncomboBox_currentIndexChanged();

    void on_PressureIncomboBox_currentIndexChanged();

    void on_TemperatureOutcomboBox_currentIndexChanged();

    void on_PressureOutcomboBox_currentIndexChanged();

    void on_resetButton_clicked();


private:
    Ui::WaterSteam *ui;
};

#endif // WATERSTEAM_H
