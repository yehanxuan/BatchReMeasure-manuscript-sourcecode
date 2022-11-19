// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
arma::uvec my_setdiff(arma::uvec& x, const arma::uvec& y){
  
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y[j]));
    x.shed_row(q1);
  }
  return x;
}

// [[Rcpp::export]]
double Update_rho_S1(mat Zc1, mat Zc2, vec Yc1, vec Yc2, double a0H, double a1H,
                  vec betaH, double sigma1H, double sigma2H, uvec Index) {
  int nc2 = Zc2.n_rows;
  vec mean1 = Zc1 * betaH;
  vec mean3 = a1H + Zc2 * betaH;
  Index = Index - 1;  // Index different from R
  
  vec Ys1 = Yc1.elem(Index);
  double W1p = sum( (Ys1 - mean1.elem(Index))%(Ys1 - mean1.elem(Index)));
  double W3p = sum( (Yc2 - mean3)% (Yc2 - mean3)); 
  double W13p = sum( (Ys1 - mean1.elem(Index))%(Yc2 - mean3) );
  
  double c3 = 1, c2 = -W13p/(nc2*sigma1H*sigma2H);
  double c1 = (W1p/(nc2*sigma1H*sigma1H) +  W3p/(nc2*sigma2H*sigma2H) - 1),
    c0 = -W13p/(nc2*sigma1H*sigma2H);
  vec P(4);
  P(0) = c3, P(1) = c2, P(2) = c1, P(3) = c0;
  cx_vec R = roots(P);
  double rhoH = real(R[2]); // real number 
  if (rhoH > 0.99)
    rhoH = 0.97;
  else if (rhoH <-0.99)
    rhoH = -0.97; 
  return rhoH;
}

// [[Rcpp::export]] 
double Update_sigma1_S1(mat Zc1, mat Zc2, vec Yc1, vec Yc2, double a0H, double a1H,
                     vec betaH, double rhoH, double sigma2H, uvec Index) {
  int nc1 = Yc1.size(), nc2 = Yc2.size();
  vec mean1 = Zc1 * betaH;
  vec mean3 = a1H + Zc2 * betaH;
  Index = Index - 1;
  uvec allid = regspace<uvec>(0, 1, nc1 - 1);
  uvec idx_c = my_setdiff(allid, Index);
  
  vec Ys1 = Yc1.elem(Index);
  double W1p = sum( (Ys1 - mean1.elem(Index))%(Ys1 - mean1.elem(Index)));
  double W1s;
  if (nc2 == nc1)
     W1s = 0;
  else 
     W1s = sum( (Yc1.elem(idx_c) - mean1.elem(idx_c)) % 
                      (Yc1.elem(idx_c) - mean1.elem(idx_c)));
  double W13p = sum( ( Yc1.elem(Index) - mean1.elem(Index) )%(Yc2 - mean3));
  
  if (rhoH > 0.99) 
    rhoH = 0.97;
  else if (rhoH < -0.99)
    rhoH = -0.97;
  double c0 = -W1p/nc1 - (1 - pow(rhoH,2))*W1s/nc1, c1 = rhoH*W13p/(nc1*sigma2H),
    c2 = (1 - pow(rhoH, 2));
  vec P = {c2, c1, c0};
  cx_vec R=roots(P);
  vec RT = real(R);
  double sigma1H = as_scalar(  RT.elem(find(RT > 0 )) );
  return sigma1H;
}

// [[Rcpp::export]]
double Update_sigma2_S1(mat Zc1, mat Zt2, mat Zc2, vec Yc1, vec Yt2, vec Yc2,
                     double a0H, double a1H,
                     vec betaH, double rhoH, double sigma1H, uvec Index) {
  int nc2 = Zc2.n_rows; 
  int nt2 = Zt2.n_rows; 
  int N2 = nc2 + nt2;
  vec mean1 = Zc1*betaH; 
  vec mean2 = a0H + a1H + Zt2 * betaH;
  vec mean3 = a1H + Zc2 * betaH;
  
  Index = Index - 1;
  double W2 = sum( (Yt2 - mean2)% (Yt2 - mean2));
  double W3p = sum((Yc2 - mean3)%(Yc2 - mean3) );
  double W13p = sum( (Yc1.elem(Index) - mean1.elem(Index)) % (Yc2 - mean3) );
  if (rhoH > 0.99)
    rhoH = 0.97;
  else if (rhoH < -0.99)
    rhoH = -0.97;
  
  double c2 = (1 - pow(rhoH, 2) ), c1 = rhoH*W13p/(N2*sigma1H),
    c0 = -W3p/N2 -(1 - pow(rhoH, 2) )*W2/N2;
  
  vec P = {c2, c1, c0};
  cx_vec R=roots(P);
  vec RT = real(R);
  double sigma2H = as_scalar(  RT.elem(find(RT > 0 )) );
  return sigma2H;
}

// [[Rcpp::export]]
vec Update_beta_S1(mat Zc1, mat Zt2, mat Zc2, vec Yc1, vec Yt2, vec Yc2, 
                   double a0H, double a1H, double rhoH, double sigma1H, double sigma2H, uvec Index) {
  int nc1 = Yc1.size();
  int nc2 = Yc2.size();
  double R = rhoH*sigma2H/sigma1H;
  Index = Index - 1;
  uvec total_vec = linspace<uvec>(0, nc1-1, nc1);
  mat Cov1 = Zc2.t() * Zc2/( pow(sigma1H,2)*(1- pow(rhoH, 2)));
  mat Cor1 = Zc2.t() * Yc1.elem(Index)/( pow(sigma1H,2)*(1- pow(rhoH, 2)));
  
  uvec allid = regspace<uvec>(0, 1, nc1 - 1);
  uvec idx_c = my_setdiff(allid, Index);
  
  if (nc2 != nc1) {
    mat Zc1cs = Zc1.rows(idx_c);
    Cov1 = Cov1 + Zc1cs.t() * Zc1cs/(sigma1H*sigma1H);
    Cor1 = Cor1 + Zc1cs.t() * Yc1.elem(idx_c)/(sigma1H*sigma1H);
  }
  
  rowvec Zc2mean = mean(Zc2, 0);
  mat Ztilde = Zc2.each_row() - (1-R)*Zc2mean;
  mat Zt2_ct = Zt2.each_row() - mean(Zt2, 0);
  
  mat  Cov2 = ( -2*rhoH*Zc2.t()*Ztilde/(sigma1H*sigma2H) + 
    Ztilde.t() * Ztilde/(sigma2H*sigma2H) )/(1 - rhoH*rhoH) + Zt2_ct.t() * Zt2_ct/(sigma2H*sigma2H);
  
  vec Ytilde = Yc2 - mean(Yc2) + R*mean(Yc1.elem(Index));
  
  vec Cor2 = ( -Zc2.t()*Ytilde - Ztilde.t()*Yc1.elem(Index) )*rhoH/(sigma1H*sigma2H*(1-rhoH*rhoH)) + 
    Ztilde.t()*Ytilde/(sigma2H*sigma2H*(1-rhoH*rhoH)) + Zt2_ct.t() * (Yt2 - mean(Yt2) )/(sigma2H*sigma2H);
  
  mat S = Cov1 + Cov2; 
  vec t = Cor1 + Cor2; 
  
  
  vec betaH = solve(S, t);
  
  //mat Zc2-Zc2mean; 
  
  //find(Index, total_vec);
 // find(totol_vec == Index);
 return betaH;
}



// [[Rcpp::export]]
vec Update_a_S1(mat Zc1, mat Zt2, mat Zc2, vec Yc1, vec Yt2, vec Yc2, 
                   double a0H, double a1H, vec betaH, double rhoH, double sigma1H, double sigma2H, uvec Index) {

  double R = rhoH*sigma2H/sigma1H;
  Index = Index - 1;
  //
  a1H = mean( Yc2 - R*Yc1.elem(Index) - (Zc2 - R*Zc1.rows(Index))*betaH );
  a0H = mean(Yt2 - Zt2*betaH) - a1H;
  vec a(2);
  a(0) = a1H; a(1) = a0H;
  return a;
}




// [[Rcpp::export]]
double Rcpp_Objective_S1(mat Zc1,mat Zt2,mat Zc2, vec Yc1, vec Yt2, vec Yc2,
                    double a0H, double a1H, vec betaH, double rhoH, double sigma1H, double sigma2H,uvec Index) {
  // uvec unsigned integer 
  Index = Index - 1;
  int nc1 = Yc1.size();
  int nc2 = Yc2.size();
  int nt2 = Yt2.size();
  uvec allid = regspace<uvec>(0, 1, nc1 - 1);
  uvec idx_c = my_setdiff(allid, Index);
  
  vec Yc1Scale = (Yc1 - Zc1*betaH)/sigma1H;
  vec Yt2Scale = (Yt2 - a0H - a1H - Zt2*betaH)/sigma2H;
  vec Yc2Scale = (Yc2 - a1H - Zc2*betaH)/sigma2H;
  
  double part1 = as_scalar( nc2*log(sigma1H) + nc2*log(sigma2H) + nc2*log(1 - pow(rhoH, 2))/2 + 
    1/(2*(1-pow(rhoH, 2))) * ( Yc1Scale.elem(Index).t()*Yc1Scale.elem(Index) - 
    2*rhoH*Yc1Scale.elem(Index).t() * Yc2Scale +  Yc2Scale.t() * Yc2Scale 
    ) + nt2*log(sigma2H) + Yt2Scale.t()*Yt2Scale/2 );
  
  double part2;
  if (nc2 == nc1) 
    part2 = 0;
  else 
    part2 = as_scalar( (nc1 - nc2)*log(sigma1H) + Yc1Scale.elem(idx_c).t() * Yc1Scale.elem(idx_c)/2 );
  
  return part1 + part2;
  
}




// /***R
// Update_a_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, sigma2H, Index)
// ***/






