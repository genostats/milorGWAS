#include <RcppEigen.h>
#include <cmath>
#include <iostream>
#ifndef GASTONLOGITMODELSCORE
#define GASTONLOGITMODELSCORE

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

using namespace Rcpp;
using namespace Eigen;


template<typename scalar_t>
scalar_t square(scalar_t x) {
  return x*x;
}


// Test du score pour la régression logistique
// On a fait le fit Y ~ X et on demande si G a un effet
//
// in:
// y1 = vecteur y - pi avec pi = expit(X beta0), calculé sous H0 (pas d'effet de G)
// G  = vecteur à tester
// W  = vecteur des pi(1-pi)
// A  = matrice telle que A'A = WX (X'WX)^{-1} X'W
// out:
// score = la stat du score
// variance = sa variance
template<typename scalar_t>
void logistic_model_score(const VECTOR<scalar_t> & y1, const VECTOR<scalar_t> & G, const VECTOR<scalar_t> & W, const MATRIX<scalar_t> & A,
                          scalar_t & score, scalar_t & variance) {
  auto sq = std::ptr_fun(square<scalar_t>);
  score = G.dot(y1);
  variance = G.unaryExpr(sq).cwiseProduct(W).sum() - (A * G).unaryExpr(sq).sum();
}

// idem quand X est réduit à l'intercept
// le calcul est plus rapide 
// Y1 et G = idem ci dessous
// w = pi*(1 - pi) (un scalaire car c'est la même valeur pour tout i)
template<typename scalar_t>
void logistic_model_score_nocovar(const VECTOR<scalar_t> & y1, const VECTOR<scalar_t> & G, scalar_t w, 
                          scalar_t & score, scalar_t & variance, bool compute_variance = true) {
  score = G.dot(y1);
  if(compute_variance) {
    auto sq = std::ptr_fun(square<scalar_t>);
    variance = w * ( G.unaryExpr(sq).sum() - square<scalar_t>(G.sum()) / (scalar_t) G.size());
  }
}

#endif
