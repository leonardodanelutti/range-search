#include <QApplication>
#include "visualization.h"

int main(int argc, char* argv[]) {
    // Qt GUI Application
    QApplication app(argc, argv);
    
    RangeSearchMainWindow window;
    window.show();
    
    return app.exec();
}
