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

private:
    Ui::WaterSteam *ui;
};

#endif // WATERSTEAM_H
