#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppEigenForward.h>
#include <RcppEigenWrap.h>
//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;
using namespace Rcpp;



// [[Rcpp::export]]
Eigen::MatrixXd inverse_cpp(Eigen::MatrixXd &objs){
  Eigen::MatrixXd obj_inverse(objs.inverse());
  return obj_inverse;
}


// [[Rcpp::export]]
Eigen::MatrixXd lag_transform(Eigen::MatrixXd objs, int t, int n, int lag,
                              bool top){
  int k(objs.cols());
  objs.resize(t,k*n);
  Eigen::MatrixXd lagM(t-lag,n*k);
  if(top == TRUE){
    lagM = objs.topRows(t-lag).array();
  }else{
    lagM = objs.bottomRows(t-lag).array();
  }
  lagM.resize((t-lag)*n,k);
  return lagM;
}

// [[Rcpp::export]]
double three_one(Eigen::VectorXd &pars1,Eigen::MatrixXd &delty0,
                 Eigen::MatrixXd &evs,Eigen::MatrixXd &omega, int &cd, int &tt,
                 int &nn,int &ny,double &varv1){

  double likely_function_value,ww,tn,p3,p2,p20,nb,sg,p1;

  Eigen::MatrixXd omega0(tt-1,tt-1),osolve(tt-1,tt-1),xn(tt-1,1),yn(tt-1,1),fits(tt-1,1),
  betas(cd-1,1),pp(1,1);

  int t;
  Eigen::VectorXd ys(ny);

  ys(0) = pars1(0);
  for(t=1;t<ny;t++){
    ys(t) = pars1(t)+ys(t-1);
  }

  double gamma = (ys.array().abs()).eval().maxCoeff();

  sg = pars1(cd-1);

  ww = varv1/sg;

  betas.col(0) = pars1.head(cd-1);

  tn = tt-1;

  nb = nn;

  p1 = tn*nb*log(6.3);

  omega0 = omega;

  if(sg<=0){
    likely_function_value = numeric_limits<double>::infinity();
  }else if(gamma > 1){
    likely_function_value = numeric_limits<double>::infinity();
  }else if(ww <= 1){
    likely_function_value = numeric_limits<double>::infinity();
  }else{

    omega0(0,0) = ww;
    omega0 = omega0*sg;
    pp(0,0) = 0;
    osolve = omega0.inverse();
    p3 = 0;

    for(t=0;t<nn;t++){
      xn = evs.middleRows(t*(tt-1),tt-1);
      yn = delty0.middleRows(t*(tt-1),tt-1);
      fits = (yn - xn*betas);
      pp = fits.adjoint()*osolve*fits;
      p3 += pp(0,0);
    }

    p20 = (1 + (tn-1)*(ww-1))*pow(sg,(tt-1));
    p2 = nb*log(p20);
    likely_function_value = 0.5*(p1+p2+p3);
    if(likely_function_value<0){
      likely_function_value = numeric_limits<double>::infinity();
    }

  }

  return likely_function_value;
 // return List::create(_["p1"] =p1,_["p2"] =p2,_["p3"] =p3,_["ww"] =ww,_["gamma"] =gamma,
 //                     _["sg"] =sg);

}

// [[Rcpp::export]]
double three_oneb(Eigen::VectorXd &pars1,Eigen::MatrixXd &delty0,
                 Eigen::MatrixXd &evs,Eigen::MatrixXd &omega, int &cd, int &tt,
                 int &nn,int &ny,double &varv1){

  double likely_function_value,ww,tn,p3,p2,p20,nb,sg,p1;

  Eigen::MatrixXd omega0(tt-1,tt-1),osolve(tt-1,tt-1),xn(tt-1,1),yn(tt-1,1),fits(tt-1,1),
  betas(cd-1,1),pp(1,1);

  int t;
  Eigen::VectorXd ys(ny);

  ys(0) = pars1(0);
  for(t=1;t<ny;t++){
    ys(t) = pars1(t);
  }

  double gamma = (ys.array().abs()).eval().maxCoeff();

  sg = pars1(cd-1);

  ww = varv1/sg;

  betas.col(0) = pars1.head(cd-1);

  tn = tt-1;

  nb = nn;

  p1 = tn*nb*log(6.3);

  omega0 = omega;

  if(sg<=0){
    likely_function_value = numeric_limits<double>::infinity();
  }else if(gamma > 1){
    likely_function_value = numeric_limits<double>::infinity();
  }else if(ww <= 1){
    likely_function_value = numeric_limits<double>::infinity();
  }else{

    omega0(0,0) = ww;
    omega0 = omega0*sg;
    pp(0,0) = 0;
    osolve = omega0.inverse();
    p3 = 0;

    for(t=0;t<nn;t++){
      xn = evs.middleRows(t*(tt-1),tt-1);
      yn = delty0.middleRows(t*(tt-1),tt-1);
      fits = (yn - xn*betas);
      pp = fits.adjoint()*osolve*fits;
      p3 += pp(0,0);
    }

    p20 = (1 + (tn-1)*(ww-1))*pow(sg,(tt-1));
    p2 = nb*log(p20);
    likely_function_value = 0.5*(p1+p2+p3);
    if(likely_function_value<0){
      likely_function_value = numeric_limits<double>::infinity();
    }

  }

  return likely_function_value;
  // return List::create(_["p1"] =p1,_["p2"] =p2,_["p3"] =p3,_["ww"] =ww,_["gamma"] =gamma,
  //                     _["sg"] =sg);

}
// [[Rcpp::export]]
double three_two(Eigen::VectorXd &pars1,Eigen::MatrixXd &delty0,
                 Eigen::MatrixXd &evs,Eigen::MatrixXd &omega, int &cd, int &tt, int &nn){
  double likely_function_value,ww,tn,p3,p2,p20,nb,sg,p1;
  Eigen::MatrixXd omega0(tt-1,tt-1),osolve(tt-1,tt-1),xn(tt-1,1),yn(tt-1,1),fits(tt-1,1),
  betas(cd-2,1),pp(1,1);
  int t;

  sg = pars1(cd-2);
  ww = pars1(cd-1);
  betas.col(0) = pars1.head(cd-2);
  tn = tt-1;
  nb = nn;
  p1 = tn*nb*log(6.3);
  omega0 = omega;
  if(sg<=0){
    likely_function_value = numeric_limits<double>::infinity();
  }else if(ww <= 1){
    likely_function_value = numeric_limits<double>::infinity();
  }else{

    omega0(0,0) = ww;
    omega0 = omega0*sg;
    pp(0,0) = 0;
    osolve = omega0.inverse();
    p3 = 0;

    for(t=0;t<nn;t++){
      xn = evs.middleRows(t*(tt-1),tt-1);
      yn = delty0.middleRows(t*(tt-1),tt-1);
      fits = (yn - xn*betas);
      pp = fits.adjoint()*osolve*fits;
      p3 += pp(0,0);
    }

    p20 = (1 + (tn-1)*(ww-1))*pow(sg,(tt-1));
    p2 = nb*log(p20);
    likely_function_value = 0.5*(p1+p2+p3);
    if(likely_function_value<0){
      likely_function_value = numeric_limits<double>::infinity();
    }

  }

  return likely_function_value;

}
