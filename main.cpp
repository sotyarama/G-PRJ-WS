#include "watersteam.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    WaterSteam w;
    w.show();

    return a.exec();
}
