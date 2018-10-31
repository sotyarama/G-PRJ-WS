#ifndef WATERSTEAM_H
#define WATERSTEAM_H

#include <QMainWindow>

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


private:
    Ui::WaterSteam *ui;
};

#endif // WATERSTEAM_H
