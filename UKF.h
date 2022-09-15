#ifndef UKF_H
#define UKF_H

/* 
	Author: GhasemHeydari 
	Email: 	ghasem.heydari2@gmail.com
*/

#include <math.h>
#include "../Libraries/CMAT/Vector.h"
#include "../Libraries/CMAT/Matrix.h"
#include "../Libraries/CMAT/MATLAB.h"
#include "../Libraries/SIM/SIMTIME.h"
#include "../Libraries/SIM/SOLVER.h"
#include "../Navigation/INS.h"
class Filter;
class Vector;
class Matrix;
class SimTime;
class Solver;
class INS;

class UKF: public StateSpaceModel
{
public:
    UKF(Vector EstimateStates, Vector Measurments, double alpha, double beta, double kapa);
    ~UKF();
    void Predict(Vector EstimateStates, Matrix Covariance, Matrix ProccessNoise, Vector InputAccelerations, Vector InputAngularrate);
    void Update(Vector EstimateStates, Matrix Covariance, Matrix MesurmentNoise, Vector InputCurrentMeasurment);
    void Step();
    Vector GetoutputStates();
    Matrix GetOutputCovariance();

private:
    Matrix CalcSigmaPoint(Vector States, Matrix S, double gamma);
    Vector CalcMeasurementEstimate(Vector Point);
    void output(SimTime &UKFtime);

private:
    /**
     * @param L: Number of states
     */
    double L;
    /**
     * @param m: Number of measurments
     */
    double m;
    /**
     * @param  Scaling factors
     */
    double lambda, c, gamma;
    /**
     * @param Wm: Weight vectors
     */
    Vector Wm, Wc;
    Matrix DiagonalWC;
    Vector xHat;
    Matrix PMinus;
    Matrix PPlus;

private:
    SimTime* UKFtime;
    Solver* rungeKutta;

private:
    /**
     * @brief Dynamic is simulation in the predict part and would need to set with your dynamic
     */
    INS* Dynamic;
};

#endif // UKF_H
