#include "pricer.h"
#include "ui_pricer.h"

#include "method.h"
#include "newton_raphson.h"

using namespace std;

Pricer::Pricer(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Pricer)
{
    ui->setupUi(this);
}

Pricer::~Pricer()
{
    delete ui;
}


void Pricer::on_tabWidget_tabBarClicked(int index)
{
    switch(index){
        //case(0): //european
            /*QMessageBox msgBox;
            msgBox.setText("The document has been modified.");
            msgBox.exec();*/
        /*case(1): // american 
        case(2): // g_asian (1)should set default value of M and n
        case(3): // a_asian (1)should set default value of M and n
        case(4): // g_basket (1)should set default value of M
        case(5): // a_basket (1)should set default value of M
        case(6): // vol*/
    }
}

void Pricer::on_button_euro_clicked()
{
    
    double T = QString(ui->expiryTimeLineEdit->text()).toDouble();
    double S = QString(ui->s0LineEdit->text()).toDouble();
    double K = QString(ui->strikePriceLineEdit->text()).toDouble();
    double r = QString(ui->nonRiskInterestRateLineEdit->text()).toDouble();
    double q = QString(ui->repoRateLineEdit->text()).toDouble();
    double sigma = QString(ui->volatilityLineEdit->text()).toDouble();
    string type = (ui->optionTypeComboBox->currentText()).toStdString();
    cout << S << "," <<K<<","<<T<<","<<r<<","<<q<<","<<sigma<<","<<type << endl;

    BlackScholes* option = new BlackScholes(S, K, r, q, T, type);
    double val = option->option_price(sigma);
    delete option;
    cout  << val << endl;
    ui->res_euro->setText(QString::number(val));
}



void Pricer::on_button_amri_clicked()
{
    double T = QString(ui->expiryTime_amri->text()).toDouble();
    double S = QString(ui->s0_amri->text()).toDouble();
    double K = QString(ui->strikePrice_amri->text()).toDouble();
    double r = QString(ui->nonRiskInterestRate_amri->text()).toDouble();
    double sigma = QString(ui->volatility_amri->text()).toDouble();
    string type = (ui->optionType_amri->currentText()).toStdString();
    int N = (ui->sizeBinomial->value());
    cout << S << "," <<K<<","<<T<<","<<r<<","<<N<<","<<sigma<<","<<type << endl;

    Binomial* amriOption = new Binomial(S, sigma, r, T, K, type, N);
    double val = amriOption->OptionPrice();
    delete amriOption;
    cout  << val << endl;
    ui->res_amri->setText(QString::number(val));
}

void Pricer::on_button_iv_clicked()
{
    double T = QString(ui->expiryTime_iv->text()).toDouble();
    double S = QString(ui->s0_iv->text()).toDouble();
    double K = QString(ui->strikePrice_iv->text()).toDouble();
    double r = QString(ui->nonRiskInterestRate_iv->text()).toDouble();
    double q = QString(ui->repoRate_iv->text()).toDouble();
    string type = (ui->optionType_iv->currentText()).toStdString();
    double C_M = QString(ui->marketPrice->text()).toDouble();
    cout << S << "," <<K<<","<<T<<","<<r<<","<<q<<","<<C_M<<","<<type << endl;
    BlackScholes bsc(S, K, r, q, T, type); 

    double init = 0.3;  // Our guess impl. vol of 30%
    double epsilon = 0.00001;

    // Calculate the implied volatility
    double sigma = newton_raphson<BlackScholes, 
                                &BlackScholes::option_price,
                                &BlackScholes::option_vega>
    (C_M, init, epsilon, bsc);

    // Output the values
    cout << "Implied Vol: " << sigma << endl;
    ui->res_iv->setText(QString::number(sigma));
}

void Pricer::on_button_ga_clicked()
{
    double T = QString(ui->expiryTime_ga->text()).toDouble();
    double S = QString(ui->s0_ga->text()).toDouble();
    double K = QString(ui->strikePrice_ga->text()).toDouble();
    double r = QString(ui->nonRiskInterestRate_ga->text()).toDouble();
    double sigma = QString(ui->volatility_ga->text()).toDouble();
    string type = (ui->optionType_ga->currentText()).toStdString();
    int n = ui->obsTimes_ga->value();
    int M = ui->M_ga->value();

    Asian* geom_asian = new Asian(S, sigma, r, T, K, n, type);
    double valByFormula = geom_asian->Geo();
    cout<<"close-form formula: " << valByFormula << endl;
    geom_asian->MonteCarlo(M);
    Conf mc = geom_asian->getMC();
    cout << "Monte Carlo: " << mc.value << ", [" << mc.confs[0] << "," << mc.confs[1] << "]" << endl;
    delete geom_asian;
    ui->res_ga->setText("Value solved by close-form formula: "+QString::number(valByFormula)+"\n"+"Value solved by Monte Carlo: "+QString::number(mc.value)+"\nConfidence intervals: ["+QString::number(mc.confs[0])+", "+QString::number(mc.confs[1])+"]");
}

void Pricer::on_button_aa_clicked()
{
    double T = QString(ui->expiryTime_aa->text()).toDouble();
    double S = QString(ui->s0_aa->text()).toDouble();
    double K = QString(ui->strikePrice_aa->text()).toDouble();
    double r = QString(ui->nonRiskInterestRate_aa->text()).toDouble();
    double sigma = QString(ui->volatility_aa->text()).toDouble();
    string type = (ui->optionType_aa->currentText()).toStdString();
    int n = ui->obsTimes_aa->value();
    int M = ui->M_aa->value();

    Arith_Asian* arith_asian = new Arith_Asian(S, sigma, r, T, K, n, type);
    arith_asian->MonteCarlo(M);
    Conf mc = arith_asian->getMC();
    cout << "Monte Carlo: " << mc.value << ", [" << mc.confs[0] << "," << mc.confs[1] << "]" << endl;
    Conf cv = arith_asian->getCV();
    cout << "Control Variable: " << cv.value << ", [" << cv.confs[0] << "," << cv.confs[1] << "]" << endl;
    delete arith_asian;
    ui->res_aa->setText("Value solved by Monte Carlo: "+QString::number(mc.value)+"\nConfidence intervals: ["+QString::number(mc.confs[0])+", "+QString::number(mc.confs[1])+"]"+"\nValue solved by Control Variable: "+QString::number(cv.value)+"\nConfidence intervals: ["+QString::number(cv.confs[0])+", "+QString::number(cv.confs[1])+"]");
}

void Pricer::on_button_gb_clicked()
{
    double T = QString(ui->expiryTime_gb->text()).toDouble();
    double S1 = QString(ui->s1_gb->text()).toDouble();
    double S2 = QString(ui->s1_gb->text()).toDouble();
    double p = QString(ui->p_gb->text()).toDouble();
    double K = QString(ui->strikePrice_gb->text()).toDouble();
    double r = QString(ui->nonRiskInterestRate_gb->text()).toDouble();
    double sigma1 = QString(ui->sigma1_gb->text()).toDouble();
    double sigma2 = QString(ui->sigma2_gb->text()).toDouble();
    string type = (ui->optionType_gb->currentText()).toStdString();
    int M = ui->M_gb->value();
    // cout<<K<<","<<r<<","<<type<<","<<M<<endl;
    int n = 2; // default value n = 2
    double S[] = {S1, S2};
    double sigma[] = {sigma1, sigma2};
    double** pp = new double*[n];
    for(int i=0; i<n; i++){
        // cout<<S[i]<<","<<sigma[i]<<endl;
        pp[i] = new double[n];
        for(int j=0; j<n; j++){ 
            if(j==i){ pp[i][j] = 1.;}
            else{ pp[i][j] = p;}
            cout<<pp[i][j]<<endl;
        }
    }

    Basket* geom_basket = new Basket(S, sigma, pp, r, T, K, n, type);
    double valByFormula = geom_basket->Geo();
    cout<<"close-form formula: " << valByFormula << endl;
    geom_basket->MonteCarlo(M);
    Conf mc = geom_basket->getMC();
    cout << "Monte Carlo: " << mc.value << ", [" << mc.confs[0] << "," << mc.confs[1] << "]" << endl;
    delete geom_basket;
    ui->res_gb->setText("Value solved by close-form formula: "+QString::number(valByFormula)+"\n"+"Value solved by Monte Carlo: "+QString::number(mc.value)+"\nConfidence intervals: ["+QString::number(mc.confs[0])+", "+QString::number(mc.confs[1])+"]");
}

void Pricer::on_button_ab_clicked()
{
    double T = QString(ui->expiryTime_ab->text()).toDouble();
    double S1 = QString(ui->s1_ab->text()).toDouble();
    double S2 = QString(ui->s2_ab->text()).toDouble();
    double p = QString(ui->p_ab->text()).toDouble();
    double K = QString(ui->strikePrice_ab->text()).toDouble();
    double r = QString(ui->nonRiskInterestRate_ab->text()).toDouble();
    double sigma1 = QString(ui->sigma1_ab->text()).toDouble();
    double sigma2 = QString(ui->sigma2_ab->text()).toDouble();
    string type = (ui->optionType_ab->currentText()).toStdString();
    int M = ui->M_ab->value();

    int n = 2;
    double S[] = {S1, S2};
    double sigma[] = {sigma1, sigma2};
    double** pp = new double*[n];
    for(int i=0; i<n; i++){
         cout<<S[i]<<","<<sigma[i]<<endl;
        pp[i] = new double[n];
        for(int j=0; j<n; j++){ 
            if(j==i){ pp[i][j] = 1.;}
            else{ pp[i][j] = p;}
            cout<<pp[i][j]<<endl;
        }
    }
    Arith_Basket* arith_basket = new Arith_Basket(S, sigma, pp, r, T, K, n, type);
    arith_basket->MonteCarlo(M);
    Conf mc = arith_basket->getMC();
    cout << "Monte Carlo: " << mc.value << ", [" << mc.confs[0] << "," << mc.confs[1] << "]" << endl;
    Conf cv = arith_basket->getCV();
    cout << "Control Variable: " << cv.value << ", [" << cv.confs[0] << "," << cv.confs[1] << "]" << endl;
    delete arith_basket;
    ui->res_ab->setText("Value solved by Monte Carlo: "+QString::number(mc.value)+"\nConfidence intervals: ["+QString::number(mc.confs[0])+", "+QString::number(mc.confs[1])+"]"+"\nValue solved by Control Variable: "+QString::number(cv.value)+"\nConfidence intervals: ["+QString::number(cv.confs[0])+", "+QString::number(cv.confs[1])+"]");
}


