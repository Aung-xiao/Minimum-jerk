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
    MatrixXd A=MatrixXd::Zero((p_num1d+1)*m, p_num1d*m);
    MatrixXd P=MatrixXd::Zero(p_num1d*m, 1);
    MatrixXd Deq=MatrixXd::Zero(p_num1d*(m+1), 3);

    for(int i=0;i<m;i++)
    {
//      double T_start=(i>0?Time(i-1):0);
      MatrixXd q=QGeneration(p_num1d,Time(i));
      MatrixXd a=AGeneration(p_num1d,d_order,Time(i));
//      MatrixXd a_next=AGeneration(p_num1d,d_order,T_next);
//      MaxtrixInsert(Q,q,(i-1)*p_num1d,(i-1)*p_num1d);
//      if(i==1)
//        MaxtrixInsert(A,a,0,0);
//      else if(i==m)
//        MaxtrixInsert(A,a_next,i*p_num1d,i*p_num1d);
//      else
//      {
//        MaxtrixInsert(A,a,i*p_num1d,(i-1)*p_num1d);
//        MaxtrixInsert(A,-a_next,i*p_num1d,i*p_num1d);
//      }
    }
//    P_=P.sparseView();
//    A_=A.sparseView();
//
//    for(int i=0;i<3;i++)
//    {
//      MatrixXd deq_single=DeqGeneration(Path.block(0,i,m,1),
//                                          Vel.block(0,i,2,1),
//                                          Acc.block(0,i,2,1),
//                                          p_num1d,m);
//      l_d=deq_single;
//      u_d=deq_single;
//      qpSolver_.setMats(P_, q_d, A_, l_d, u_d);
//      qpSolver_.solve();
//      int ret = qpSolver_.getStatus();
//      if (ret != 1)
//        ROS_ERROR("fail to solve QP!");
//      Eigen::VectorXd sol = qpSolver_.getPrimalSol();
//      MaxtrixInsert(PolyCoeff,sol,0,i);
//    }

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
  cout<<"q:"<<endl;
  cout<<Q<<endl;
  return Q;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::AGeneration(int p_num1d, int d_order, double T_end) {
  MatrixXd A= MatrixXd::Zero(p_num1d,  p_num1d);
  for(int i=d_order;i<p_num1d;i++)
  {
    A(i,i)=fac(i)/fac(i-d_order)*pow(T_end,i-d_order);
  }
  cout<<"a:"<<endl;
  cout<<A<<endl;
  return A;
}

int TrajectoryGeneratorWaypoint::fac(int x) {
    int f;
    if(x==0 || x==1)
      f=1;
    else
      f=fac(x-1)*x;
    return f;
}
void TrajectoryGeneratorWaypoint::MaxtrixInsert(Eigen::MatrixXd& M,Eigen::MatrixXd m,int row,int col)
{
  for(int i=row+m.rows();i>=m.rows();i--)
  {
    for(int l=col+m.cols();l>=m.cols();l--)
    {
      M(i,l)=m(i-m.rows(),l-m.cols());
    }
  }
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::DeqGeneration(Eigen::MatrixXd Path,Eigen::MatrixXd Vel,Eigen::MatrixXd Acc,int p_num1d,int m)
{
  MatrixXd Deq=MatrixXd::Zero(p_num1d*(m+1),1);
  int l=0;
  for(int i=0;i<p_num1d*m;i+=p_num1d)
  {
    Deq(i+p_num1d/2+1,1)=Path(l);
    l++;
  }
  Deq(0+p_num1d/2+2,1)=Vel(0);
  Deq(0+p_num1d/2+3)=Acc(0);
  Deq(p_num1d*(m-1)+p_num1d/2+2,1)=Vel(1);
  Deq(p_num1d*(m-1)+p_num1d/2+3)=Acc(1);
  return Deq;
}