#ifndef _TRAJECTORY_GENERATOR_WAYPOINT_H_
#define _TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <vector>
#include "iosqp.hpp"

class TrajectoryGeneratorWaypoint {
    private:
		double _qp_cost;
		Eigen::MatrixXd _Q;
		Eigen::VectorXd _Px, _Py, _Pz;
    public:
        TrajectoryGeneratorWaypoint();

        ~TrajectoryGeneratorWaypoint();

        Eigen::MatrixXd PolyQPGeneration(
            const int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time);
        int fac(int x);
        Eigen::MatrixXd QGeneration(int p_num1d,double T);
        Eigen::MatrixXd AConGeneration(int p_num1d,int d_order,double T_end);
        Eigen::MatrixXd AGeneration(int m, int p_num1d,int d_order,double final_T);
        Eigen::MatrixXd DeqGeneration(Eigen::MatrixXd Path,int d_order,int m);
        void matrixFill(Eigen::MatrixXd& M,Eigen::MatrixXd m,int row,int col);
        int Factorial(int x);
        Eigen::SparseMatrix<double> P_, A_;
        Eigen::VectorXd q_d,l_d,u_d;
        osqp::IOSQP qpSolver_;
};
        

#endif
