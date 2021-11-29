// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]


#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd CppMatMult(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2){
  Eigen::MatrixXd m3 = m1 * m2;
  return(m3);
}

// [[Rcpp::export]]
Eigen::MatrixXd CppMatWiseMul(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2){
  Eigen::MatrixXd m3 =  m1.cwiseProduct(m2);
  return(m3);
}

// [[Rcpp::export]]
Eigen::MatrixXd CppMatWiseDiv(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2){
  Eigen::MatrixXd m3 = m1.cwiseQuotient(m2);
  return(m3);
}

// [[Rcpp::export]]
Eigen::MatrixXd CppMatAddScalar(const Eigen::MatrixXd& m, double eps=1e-8){
  Eigen::MatrixXd mx = m;
  mx.array() += eps;
  return(mx);
}


// [[Rcpp::export]]
Eigen::MatrixXd NMF(const Eigen::MatrixXd& W, Eigen::MatrixXd Q, Eigen::MatrixXd S,  int n, int nclass, int maxiter=200, bool display_progress=true){
  double eps = 1e-8;
  Progress p(maxiter, display_progress);
  for (int i=0; i<maxiter; ++i){
      p.increment();
      // Update Q
      Eigen::MatrixXd WQS = CppMatMult(W, CppMatMult(Q, S));
      Q = CppMatWiseMul(Q, (CppMatWiseDiv(WQS, CppMatAddScalar(CppMatMult(Q, CppMatMult(Q.transpose(), WQS)), eps))).cwiseSqrt());
      // Update S
      Eigen::MatrixXd QTQ = CppMatMult(Q.transpose(), Q);
      Eigen::MatrixXd WQ = CppMatMult(W, Q);
      Eigen::MatrixXd QTWQ = CppMatMult(Q.transpose(), WQ);
      S =CppMatWiseMul(S, (CppMatWiseDiv(QTWQ, CppMatAddScalar(CppMatMult(QTQ, CppMatMult(S, QTQ)), eps))).cwiseSqrt());
  }
  return CppMatMult(Q,(S).cwiseSqrt());
}
