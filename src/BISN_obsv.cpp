#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/round.hpp>
#include "Lentz_Algorithm.hpp"
#define POSITIVE_EPS 0.0001

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List BISN_obsv(arma::mat XDat, double eta, int maxIter, double tol, double r, int s) {
  uword p = XDat.n_cols,n = XDat.n_rows;
  mat nS = XDat.t()*XDat;
  
  /* Initialization */
  double timeBegin = clock();
  vec dnS = nS.diag();
  uword pe = p*(p-1)/2;
  
  uvec idl(pe);
  uvec idu(pe);
  uvec idr(pe);
  uvec idc(pe);
  uword k = 0, i, j, jp, kappa;
  for (j = 0,jp = 0;j<p;j++,jp+=p)
  {
    for (i=j+1;i<p;i++)
    {
      idl(k) = jp+i;
      idu(k) = i*p+j;
      idr(k) = i;
      idc(k) = j;
      k++;
    }
    
  }
  
  double psr = (double)p/s;
  uvec idp(p), id0(s), idd;
  vec nd2pp = (double)n/2+linspace<vec>(p,1,p);
  
  mat ML(p,p,fill::eye);
  mat ML2(p,p,fill::eye);
  mat ML2pVL(p,p,fill::eye);
  mat VL(p,p,fill::zeros);
  mat LAMBDA(p,p,fill::zeros);
  
  vec h(pe,fill::zeros),zeta(pe),mL(pe,fill::zeros);
  zeta.fill(10.0);
  
  vec alpha(p,fill::ones);
  vec beta(p,fill::ones);
  
  double a = (double)pe/2, b = a/50;
  
  vec mLold = h/zeta;
  vec lambdaold(pe,fill::ones);
  
  mat VL2(p,p,fill::zeros);
  for (i=1;i<p;i++)
  {
    VL2(span(i,p-1),i) += i;
    VL2(i,span(i,p-1)) += i;
  }
  VL2 *= 0.01;
  
  vec d(pe,1);
  d.fill(0.5);
  
  vec mD, mD2, mD2pvD, vD, vL, lambda(pe), gmL, gvL, gmD, gvD, c5(p), c6, mLnew, lambdanew, alphatmp, betatmp, dtmp;//mL,
  mat c1, c2, c3, c4, LDL, c2pc3, K_tmp1, K_tmp2;
  double omega, difmL, diflambda, difmax;
  vec d1h(pe,fill::zeros), d1zeta(pe,fill::zeros), d1alpha(p,fill::zeros), d1beta(p,fill::zeros), d1d(pe,fill::zeros);
  double d2h = 0, d2zeta = 0, d2alpha = 0, d2beta = 0, d2d = 0, d1b = 0, d2b = 0;
  double tau = 6e2; //, tauzeta = tauh, taualpha = tauh, taubeta = tauh, taub = tauh, taud = tauh;
  double rho, gb, btmp;//1/eta;rho_ub, c7=(5*eta>0.25)?0.25:5*eta, 
  vec gh, gzeta, galpha, gbeta, gd;
  //arma_rng::set_seed(0);
  
  /* KL proximal variational inference */
  printf("Start Running BISN ...\n");
  for (kappa=1;kappa<=maxIter;kappa++)
  {
    mD = alpha/beta;
    mD2 = square(mD);
    vD = mD/beta;
    mD2pvD = mD2+vD;
    
    vL = 1/zeta;
    VL.elem(idl) = vL;
    mL =  h%vL;
    ML.elem(idl) = mL;
    ML2.elem(idl) = square(mL);
    ML2pVL.elem(idl) = ML2.elem(idl)+vL;
    
    omega = a/b;
    for (i=0;i<pe;i++) lambda(i)=(d(i)>10)?Lentz_Algorithm(d(i)):-boost::math::expint(-d(i));
    lambda.elem(find(d<=10)) %= exp(d.elem(find(d<=10)));
    lambda = 1/(d%lambda)-1;
    //         lambda.elem(find(d<0)).fill(1e6);
    //         lambda.elem(find(lambda>1e6)).fill(1e6);
    LAMBDA.elem(idl) = omega*lambda;
    LAMBDA.elem(idu) = LAMBDA.elem(idl);
    
    
    if (kappa==1)
    {
      c1 = nS;
      c2 = -VL;
      c2.each_row() += sum(VL);
      c2 *= LAMBDA(1); 
      //c2 = LAMBDA(1)*(repmat(sum(VL),p,1)-VL);
      c3 = LAMBDA;
      c4 = (VL+VL2)*mD2pvD(0);
      VL2.clear();
    }
    else
    {
      K_tmp1 = ML.rows(id0);
      K_tmp1.each_row() %= mD.t();
      LDL = K_tmp1*ML.t();
      c1 *= r;
      c1.rows(id0) += psr*(nS.rows(id0)+LDL%LAMBDA.rows(id0))*ML;
      
      c2 *= r;
      c2.rows(id0) += psr*LAMBDA.rows(id0)*VL;
      
      c3 *= r;
      c3.rows(id0) += psr*LAMBDA.rows(id0)*ML2;
      
      K_tmp1 = ML2pVL.rows(id0);
      K_tmp1.each_row() %= mD2pvD.t();
      K_tmp2 = ML2.rows(id0);
      K_tmp2.each_row() %= mD2.t();
      c4 *=r;
      c4.rows(id0) += psr*(K_tmp1*ML2pVL.t()-K_tmp2*ML2.t()+square(LDL));
    }
    
    c2pc3 = c2+c3;
    gmL = -c1.elem(idl)%mD.elem(idc)-ML.elem(idl)%(mD2pvD.elem(idc)%c2.elem(idl)+vD.elem(idc)%c3.elem(idl));
    gvL = dnS.elem(idr)%mD.elem(idc)+c2pc3.elem(idl)%mD2pvD.elem(idc);
    gh = gmL+mL%gvL - h;
    gvL.elem(find(gvL<0)).zeros();
    gzeta = gvL - zeta;
    d1h = (1-1/tau)*d1h+gh/tau;
    d2h = (1-1/tau)*d2h+mean(square(gh))/tau;
    d1zeta = (1-1/tau)*d1zeta + gzeta/tau;
    d2zeta = (1-1/tau)*d2zeta + mean(square(gzeta))/tau;
    
    gmD = (VL.t()*dnS+trans(sum(ML%c1))+trans(sum(VL%(c2pc3+c3)))%mD)/2;
    gvD = trans(sum(ML2pVL%c2pc3))/4;
    for (i=0;i<p;i++) c5(i) = boost::math::trigamma(alpha(i));
    c6 = mD/(alpha%c5-1);
    alphatmp = nd2pp+c6/beta%gvD;
    alphatmp.elem(find(alphatmp<0)).zeros();
    betatmp = gmD+(1/beta+c5%c6)%gvD;
    betatmp.elem(find(betatmp<0)).zeros();
    galpha = alphatmp - alpha;
    gbeta = betatmp - beta;
    d1alpha = (1-1/tau)*d1alpha + galpha/tau;
    d2alpha = (1-1/tau)*d2alpha + mean(square(galpha))/tau;
    d1beta = (1-1/tau)*d1beta + gbeta/tau;
    d2beta = (1-1/tau)*d2beta + mean(square(gbeta))/tau;
    
    dtmp = omega/2*c4.elem(idl);
    dtmp.elem(find(dtmp<0)).zeros();
    gd = dtmp - d;
    d1d = (1-1/tau)*d1d + gd/tau;
    d2d = (1-1/tau)*d2d + mean(square(gd))/tau;
    
    btmp = sum(lambda%c4.elem(idl))/2;
    if (btmp<0) btmp = 0;
    gb = btmp - b;
    d1b = (1-1/tau)*d1b + gb/tau;
    d2b = (1-1/tau)*d2b + pow(gb,2)/tau;
    
    rho = (mean(square(d1h))+mean(square(d1zeta))+p/pe*(mean(square(d1alpha))+mean(square(d1beta)))+mean(square(d1d))+pow(d1b,2)/pe)/(d2h+d2zeta+p/pe*(d2alpha+d2beta)+d2d+d2b/pe);
    if (rho>eta) rho = eta;
    
    
    
    
    tau = (1-rho)*tau + 1;
    h += rho*gh;
    zeta += rho*gzeta;
    alpha += rho*galpha;
    beta += rho*gbeta;
    d += rho*gd;
    b += rho*gb;
    
    
    if (kappa%100 == 0)
    {
      mLnew = h/zeta;
      lambdanew = lambda;
      difmL = sqrt(mean(square(mLnew-mLold))/mean(square(mLold)));
      diflambda = max(abs(lambdanew-lambdaold));
      difmax = max(abs(mLnew-mLold));
      printf("#no. of iterations = %d, difmL = %f, difmax = %f, diflambda = %f\n",kappa,difmL,difmax,diflambda);
      if (difmL<tol)
        break;
      else
      {
        mLold = mLnew;
        lambdaold = lambdanew;
      }
    }
    
    idp = linspace<uvec>(0, p-1, p);
    for (i=0;i<s;i++)
    {
      k = randi<uword>(distr_param(0,p-1-i));
      id0(i) = idp(k);
      if (k!=p-1-i) idp.subvec(k,p-2-i) = idp.subvec(k+1,p-1-i);
    }
    
    
    K_tmp1 = ML.rows(id0);
    K_tmp1.each_row() %= mD.t();
    LDL = K_tmp1*ML.t();
    c1.rows(id0) -= psr*(nS.rows(id0)+LDL%LAMBDA.rows(id0))*ML;
    c2.rows(id0) -= psr*LAMBDA.rows(id0)*VL;
    c3.rows(id0) -= psr*LAMBDA.rows(id0)*ML2;
    
    K_tmp1 = ML2pVL.rows(id0);
    K_tmp1.each_row() %= mD2pvD.t();
    K_tmp2 = ML2.rows(id0);
    K_tmp2.each_row() %= mD2.t();
    c4.rows(id0) -= psr*(K_tmp1*ML2pVL.t()-K_tmp2*ML2.t()+square(LDL));
  }
  
  
  
  
  ML.elem(idl) = h/zeta;
  VL.elem(idl) = 1/zeta;
  mD = alpha/beta;
  vD = mD/beta;
  for (i=0;i<pe;i++) lambda(i)=(d(i)>10)?Lentz_Algorithm(d(i)):-boost::math::expint(-d(i));
  lambda.elem(find(d<=10)) %= exp(d.elem(find(d<=10)));
  lambda = 1/(d%lambda)-1;
  double ElapsedTime = (clock()-timeBegin)/CLOCKS_PER_SEC;
  if (kappa < maxIter && difmL<tol) printf("BISN converges, elapsed time is %f seconds.\n",ElapsedTime);
  else printf("BISN reaches the maximum number of iterations, elapsed time is %f seconds.\n",ElapsedTime);
  
  
  
  return List::create(Named("ML") = ML, Named("VL") = VL, Named("mD") = mD, Named("vD") = vD,
                            Named("omega") = omega, Named("lambda") = lambda);
  
}
