#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
//RNGScope scope;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double coconetA1(Rcpp::NumericVector x, arma::mat A, arma::colvec y)
{
	int n = y.size();
	arma::colvec X = arma::colvec(n, fill::ones);
	arma::vec y_trait = y - X*x(2); 
	arma::mat y2_mat = zeros<arma::mat>(n, n);
	y2_mat.each_col() += square(y_trait); 
	arma::mat yi_yj_square = y2_mat + y2_mat.t();
	arma::mat y_trait_mat = zeros<arma::mat>(n, n);
	y_trait_mat.diag() = y_trait;
	sp_mat y_trait_mat_sp = arma::sp_mat(y_trait_mat);
	
	arma::mat Eye(n,n,fill::eye);
	arma::mat B = x(0)*Eye + x(1)*A ;
	arma::mat Bii = zeros<arma::mat>(n, n);
	Bii.each_col() += B.diag();
	arma::mat detB = square(Bii) - square(B);
	
	arma::mat Bdet = detB;
	Bdet.diag() = arma::colvec(n, fill::ones);
	double logdetB = -0.5*arma::accu(arma::vectorise(log(Bdet)));
	
	arma::mat yb = yi_yj_square%Bii-2*(y_trait_mat_sp * B) *  y_trait_mat_sp;
	double logP = 0.5/(n-1)*(logdetB -0.5*arma::accu(arma::vectorise(1.0/Bdet) % arma::vectorise( 2*yb))-0.5*n*log(2*atan(1.0)*4));
  	return -logP;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double coconetA2(Rcpp::NumericVector x, arma::mat A, arma::mat A2, arma::colvec y)
{
	int n = y.size();
	arma::colvec X = arma::colvec(n, fill::ones);
	vec y_trait = y - X*x(3); 
	mat y2_mat = zeros<arma::mat>(n, n);
	y2_mat.each_col() += arma::square(y_trait); 
	mat yi_yj_square = y2_mat + y2_mat.t();
	mat y_trait_mat = zeros<arma::mat>(n, n);
	y_trait_mat.diag() = y_trait;
	sp_mat y_trait_mat_sp = arma::sp_mat(y_trait_mat);
	
	mat Eye(n,n,fill::eye);
	mat B = x(0)*Eye + x(1)*A +x(2) * A2;
	mat Bii = zeros<mat>(n, n);
	Bii.each_col() += B.diag();
	mat detB = square(Bii) - square(B);
	
	mat Bdet = detB;
	Bdet.diag() = arma::colvec(n, fill::ones);
	double logdetB = -0.5*accu(vectorise(log(Bdet)));
	
	mat yb = yi_yj_square%Bii-2*(y_trait_mat_sp * B) *  y_trait_mat_sp;
	double logP = 0.5/(n-1)*(logdetB -0.5*arma::accu(arma::vectorise(1.0/Bdet) % arma::vectorise( 2*yb))-0.5*n*log(2*atan(1.0)*4));
  	return -logP;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double coconetA3(Rcpp::NumericVector x, arma::mat A, arma::mat A2, arma::mat A3, arma::colvec y)
{
	int n = y.size();
	arma::colvec X = arma::colvec(n, fill::ones);
	arma::vec y_trait = y - X*x(4); 
	arma::mat y2_mat = zeros<arma::mat>(n, n);
	y2_mat.each_col() += arma::square(y_trait); 
	arma::mat yi_yj_square = y2_mat + y2_mat.t();
	arma::mat y_trait_mat = zeros<arma::mat>(n, n);
	y_trait_mat.diag() = y_trait;
	arma::sp_mat y_trait_mat_sp = arma::sp_mat(y_trait_mat);
	
	arma::mat Eye(n,n,fill::eye);
	arma::mat B = x(0)*Eye + x(1)*A +x(2) * A2 + x(3) * A3 ;
	arma::mat Bii = zeros<arma::mat>(n, n);
	Bii.each_col() += B.diag();
	arma::mat detB = arma::square(Bii) - arma::square(B);
	
	arma::mat Bdet = detB;
	Bdet.diag() = arma::colvec(n, fill::ones);
	double logdetB = -0.5*arma::accu(arma::vectorise(log(Bdet)));
	
	arma::mat yb = yi_yj_square%Bii-2*(y_trait_mat_sp * B) *  y_trait_mat_sp;
	double logP = 0.5/(n-1)*(logdetB -0.5*arma::accu(arma::vectorise(1.0/Bdet) % arma::vectorise( 2*yb))-0.5*n*log(2*atan(1.0)*4));
  	return -logP;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double coconetA4(Rcpp::NumericVector x, arma::mat A, arma::mat A2,arma::mat A3, arma::mat A4, arma::colvec y)
{
	int n = y.size();
	arma::colvec X = arma::colvec(n, fill::ones);
	arma::vec y_trait = y - X*x(3); 
	arma::mat y2_mat = zeros<arma::mat>(n, n);
	y2_mat.each_col() += arma::square(y_trait); 
	arma::mat yi_yj_square = y2_mat + y2_mat.t();
	arma::mat y_trait_mat = zeros<arma::mat>(n, n);
	y_trait_mat.diag() = y_trait;
	arma::sp_mat y_trait_mat_sp = arma::sp_mat(y_trait_mat);
	
	arma::mat Eye(n,n,fill::eye);
	arma::mat B = x(0)*Eye + x(1)*A +x(2) * A2 + x(3) * A3+x(4) * A4;
	arma::mat Bii = zeros<arma::mat>(n, n);
	Bii.each_col() += B.diag();
	arma::mat detB = square(Bii) - square(B);
	
	arma::mat Bdet = detB;
	Bdet.diag() = arma::colvec(n, fill::ones);
	double logdetB = -0.5*arma::accu(arma::vectorise(log(Bdet)));
	
	arma::mat yb = yi_yj_square%Bii-2*(y_trait_mat_sp * B) *  y_trait_mat_sp;
	double logP = 0.5/(n-1)*(logdetB -0.5*arma::accu(arma::vectorise(1.0/Bdet) % arma::vectorise( 2*yb))-0.5*n*log(2*atan(1.0)*4));
  	return -logP;
}
