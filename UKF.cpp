#include "UKF.h"

// Author: 	GhasemHeydari
// Email: 	ghasem.heydari2@gmail.com

UKF::UKF(Vector EstimateStates, Vector Measurments, double alpha, double beta, double kapa)
{
    L = EstimateStates.length();
    m = Measurments.length();
    lambda = ( pow(alpha,2) * (L + kapa)) - L;
    c = L + lambda;
    Wm = Vector(2*L+1, 0.5/c);
    Wm(0) = lambda/c;
    Wc = Wm;
    Wc(0) = Wc(0) + (1 - pow(alpha,2) + beta);
    gamma = sqrt(c);

    DiagonalWC = Matrix(1 + (2 * L), 1 + (2 * L));
    for(int i = 0; i < 1 + (2 * L); i++){
        DiagonalWC(i, i) = Wc(i);
    }
    rungeKutta = Solver::create(ODE4_RungeKutta, this);
    xHat = Vector(L);
    PMinus = Matrix(L, L);
    PPlus = Matrix(L, L);
}

void UKF::Predict(Vector X, Matrix Covariance, Matrix ProccessNoise, Vector InputAccelerations, Vector InputAngularrate){

    Matrix SigmaPoints = Matrix(L, 2*L+1);
    SigmaPoints = CalcSigmaPoint(X, Covariance, gamma);

    Dynamic->InputAccelerations = InputAccelerations;
    Dynamic->InputAngularrate = InputAngularrate;

    Matrix xHatMinus = Matrix(L, 2*L+1);
    Matrix PointsOfDynamic = Matrix(L, 2*L+1);
    Matrix Y1 = Matrix(L, 2*L+1);
    Vector SigmaFromDynamic = Vector(L);

    for(int i = 0; i < SigmaPoints.cols(); i++){
        Vector x = SigmaPoints.getColumn(i);
        Dynamic->setContinuousStatesValueVector(x);
        rungeKutta->step(*(this->UKFtime));
        SigmaFromDynamic = Dynamic->getContinuousStatesValueVector();
        PointsOfDynamic.setColumn(i, SigmaFromDynamic);
        xHat += Wm(i) * SigmaFromDynamic;
        xHatMinus.setColumn(i,xHat);
    }

    Y1 = PointsOfDynamic - xHatMinus;
    PMinus = (Y1 * DiagonalWC * (Y1.transpose())) + ProccessNoise;
}

void UKF::Update(Vector EstimateStates, Matrix Covariance, Matrix MeasurmentNoise, Vector InputCurrentMeasurment){

    Matrix SigmaPoints = Matrix(L, 2*L+1);
    SigmaPoints = CalcSigmaPoint(EstimateStates, Covariance, gamma);

    Vector ZHat = Vector(m);
    Vector SigmaFromMeasurement = Vector(m);
    Matrix Zmuserment = Matrix(m, 2*L+1);
    Matrix Y2 = Matrix(m, 2*L+1);
    Matrix ZHatA = Matrix(m, 2*L+1);
    Matrix Pzz = Matrix(m, m);
    Matrix Y3 = Matrix(L, 2*L+1);;
    Matrix Pxz = Matrix(L, m);
    Matrix k = Matrix(L, m);
    Matrix xHatMinus = Matrix(L, 2*L+1);

    for(int i = 0; i < SigmaPoints.cols(); i++)
    {
        Vector x = SigmaPoints.getColumn(i);
        SigmaFromMeasurement = CalcMeasurementEstimate(x);
        Zmuserment.setColumn(i, SigmaFromMeasurement);
        ZHat += Wm(i) * SigmaFromMeasurement;
        ZHatA.setColumn(i,ZHat);
    }

    Y2 = Zmuserment - ZHatA;
    Pzz = Y2 * (DiagonalWC * Y2.transpose()) + MeasurmentNoise ;
    Y3 = SigmaPoints - xHatMinus;
    Pxz = Y3 * (DiagonalWC * Y2.transpose());
    k =  Pxz * Pzz.inverse() ;

    xHat = xHat + (k * (InputCurrentMeasurment - ZHat));
    PPlus = PMinus - (k.transpose() * Pzz * k);
}

void UKF::Step(){

}

void UKF::output(SimTime &UKFtime){
    Dynamic->output(UKFtime);
}

Matrix UKF::CalcSigmaPoint(Vector States, Matrix SCovariance, double gamma){

    int n  = States.length();
    Matrix cholesky = chol(SCovariance);
    Matrix A = gamma * cholesky.transpose();
    Matrix Y = zeros(n ,n);
    for (int i = 0; i < n; i++){
        Y.setColumn(i, States );
    }
    Matrix YplusA = Y + A;
    Matrix YminusA = Y - A;

    Matrix chi = zeros(n, ((2 * n) + 1));
    chi.setColumn(0, States);
    chi.setSubMatrix(0, 1, YplusA);
    chi.setSubMatrix(0, n + 1, YminusA);
    return chi;
}

Vector UKF::CalcMeasurementEstimate(Vector Point){

   return Point;
}

UKF::~UKF(){

}
