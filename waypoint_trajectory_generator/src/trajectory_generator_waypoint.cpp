#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

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

    MatrixXd Q=MatrixXd::Zero(p_num1d*m, p_num1d*m);
    MatrixXd A=AGeneration(m,p_num1d,d_order,Time(m-1));
    MatrixXd P=MatrixXd::Zero(p_num1d*m, 1);

    for(int i=0;i<m;i++)
    {
      MatrixXd q=QGeneration(p_num1d,Time(i));
      matrixFill(Q,q,i*p_num1d,i*p_num1d);
      if(i!=m-1)
      {
        MatrixXd a_con=AConGeneration(p_num1d,d_order,Time(i));
        matrixFill(A,a_con,(m+1)+i*d_order,i*p_num1d);
      }
    }

//    P_=P.sparseView();
//    A_=A.sparseView();
//
    for(int i=0;i<3;i++)
    {
      MatrixXd deq_s=DeqGeneration(Path.col(i),d_order,m);
//      l_d=deq_single;
//      u_d=deq_single;
//      qpSolver_.setMats(P_, q_d, A_, l_d, u_d);
//      qpSolver_.solve();
//      int ret = qpSolver_.getStatus();
//      if (ret != 1)
//        ROS_ERROR("fail to solve QP!");
//      Eigen::VectorXd sol = qpSolver_.getPrimalSol();
//      MaxtrixInsert(PolyCoeff,sol,0,i);
    }

    return PolyCoeff;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::QGeneration(int p_num1d,double T){
  MatrixXd Q= MatrixXd::Zero(p_num1d,  p_num1d);
  for(int i=3;i<p_num1d;i++)
  {
    for(int l=3;l<p_num1d;l++)
    {
      Q(i,l)=(i*(i-1)*(i-2)*l*(l-1)*(l-2))/(i+l-5)*pow(T,i+l-5);
    }
  }
//  cout<<"q:"<<Q<<endl;
  return Q;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::AConGeneration(int p_num1d, int d_order, double T_end) {
  MatrixXd A_con= MatrixXd::Zero(d_order,  2*p_num1d);
  for(int k=0;k<d_order;k++)
  {
    for(int i=k;i<p_num1d;i++)
    {
      A_con(k,i)=fac(i)/fac(i-k)*pow(T_end,i-k);
      A_con(k,i+p_num1d)=(-1)*fac(i)/fac(i-k)*pow(T_end,i-k);
    }
  }
//  cout<<"a_con:"<<A_con<<endl;
  return A_con;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::AGeneration(int m, int p_num1d, int d_order,double final_T) {
  MatrixXd A=MatrixXd::Zero((m+1)+(m-1)*d_order, p_num1d*m);
  for(int i=0;i<m;i++)
    A(i,i*p_num1d)=1;
  for(int i=0;i<p_num1d;i++)
    A(m,i+p_num1d*(m-1))=pow(final_T,i);
//  cout<<"A_Der:"<<A<<endl;
  return A;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::DeqGeneration(Eigen::MatrixXd Path,int d_order,int m)
{
  MatrixXd Deq=MatrixXd::Zero((m+1)+(m-1)*d_order,1);
  for(int i=0;i<m+1;i++)
    Deq(i,0)=Path(i);
//  cout<<"Deq:"<<Deq<<endl;
  return Deq;
}

void TrajectoryGeneratorWaypoint::matrixFill(Eigen::MatrixXd& M,Eigen::MatrixXd m,int row,int col) //非常不优雅，要改！！！
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

