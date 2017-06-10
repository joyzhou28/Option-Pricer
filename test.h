#ifndef PRICER_H
#define PRICER_H

#include <QMainWindow>
#include <QPushButton>
#include <QMessageBox>
#include <QTabWidget> 

namespace Ui {
class Pricer;
}

class Pricer : public QMainWindow
{
    Q_OBJECT

public:
    explicit Pricer(QWidget *parent = 0);
    ~Pricer();

public slots:
    void on_tabWidget_tabBarClicked(int index);
    void on_button_euro_clicked();
    void on_button_amri_clicked();
    void on_button_iv_clicked();
    void on_button_ga_clicked();
    void on_button_aa_clicked();
    void on_button_gb_clicked();
    void on_button_ab_clicked();

private:
    Ui::Pricer *ui;
};

#endif // PRICER_H
