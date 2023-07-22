#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/console.h>
#include <iostream>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments；
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients；
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    Eigen::SparseMatrix<double> Q(p_num1d*m, p_num1d*m);
    Eigen::SparseMatrix<double> A=AGeneration(m,p_num1d,d_order,Time(m-1));

    for(int i=0;i<m;i++)
    {
      Eigen::SparseMatrix<double> q=QGeneration(p_num1d,Time(i));
      spareMatrixFill(Q,q,i*p_num1d,i*p_num1d);
      if(i!=m-1)
      {
        Eigen::SparseMatrix<double> a_con=AConGeneration(p_num1d,d_order,Time(i));
        spareMatrixFill(A,a_con,(m+1)+i*d_order,i*p_num1d);
      }
    }
    VectorXd q_d(p_num1d*m);
    q_d.setZero();
    int axis=3;
    for(int l=0;l<axis;l++)
    {
      VectorXd deq_s=DeqGeneration(Path.col(l),d_order,m);
      qpSolver_.setMats(Q, q_d, A, deq_s, deq_s);
      qpSolver_.solve();
      int ret = qpSolver_.getStatus();
      if (ret != 1)
        ROS_ERROR("fail to solve QP!");
      Eigen::VectorXd sol = qpSolver_.getPrimalSol();
      for(int i=0;i<m;i++)
        matrixFill(PolyCoeff,sol.block(i*p_num1d,0,p_num1d,1).transpose(),i,l*p_num1d);
    }
    return PolyCoeff;
}

Eigen::SparseMatrix<double> TrajectoryGeneratorWaypoint::QGeneration(int p_num1d,double T){
  Eigen::SparseMatrix<double> Q(p_num1d,  p_num1d);
  for(int i=3;i<p_num1d;i++)
  {
    for(int l=3;l<p_num1d;l++)
    {
      Q.insert(i,l)=(i*(i-1)*(i-2)*l*(l-1)*(l-2))/(i+l-5)*pow(T,i+l-5);
    }
  }
  return Q;
}

Eigen::SparseMatrix<double> TrajectoryGeneratorWaypoint::AConGeneration(int p_num1d, int d_order, double T_end) {
  Eigen::SparseMatrix<double> A_con(d_order,  2*p_num1d);
  for(int k=0;k<d_order;k++)
  {
    for(int i=k;i<p_num1d;i++)
    {
      A_con.insert(k,i)=fac(i)/fac(i-k)*pow(T_end,i-k);
      A_con.insert(k,i+p_num1d)=(-1)*fac(i)/fac(i-k)*pow(0,i-k);
    }
  }
  return A_con;
}

Eigen::SparseMatrix<double> TrajectoryGeneratorWaypoint::AGeneration(int m, int p_num1d, int d_order,double final_T) {
  Eigen::SparseMatrix<double> A((m+1)+(m-1)*d_order, p_num1d*m);
  for(int i=0;i<m;i++)
    A.insert(i,i*p_num1d)=1;
  for(int i=0;i<p_num1d;i++)
    A.insert(m,i+p_num1d*(m-1))=pow(final_T,i);
  return A;
}

Eigen::VectorXd TrajectoryGeneratorWaypoint::DeqGeneration(Eigen::MatrixXd Path,int d_order,int m)
{
  Eigen::VectorXd Deq=MatrixXd::Zero((m+1)+(m-1)*d_order,1);
  for(int i=0;i<m+1;i++)
    Deq(i,0)=Path(i);
  return Deq;
}

void TrajectoryGeneratorWaypoint::spareMatrixFill(Eigen::SparseMatrix<double>& M,Eigen::MatrixXd m,int row,int col)
{
  for(int i=row;i<row+m.rows();i++)
  {
    for(int l=col;l<col+m.cols();l++)
    {
      M.insert(i,l)=m(i-row,l-col);
    }
  }
}

void TrajectoryGeneratorWaypoint::matrixFill(Eigen::MatrixXd& M,Eigen::MatrixXd m,int row,int col)
{
  for(int i=row;i<row+m.rows();i++)
  {
    for(int l=col;l<col+m.cols();l++)
    {
      M(i,l)=m(i-row,l-col);
    }
  }
}

int TrajectoryGeneratorWaypoint::fac(int x) {
  int f;
  if(x==0 || x==1)
    f=1;
  else
    f=fac(x-1)*x;
  return f;
}

