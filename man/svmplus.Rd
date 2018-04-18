\docType{package}
\name{svmplus}
\alias{svmplus}
\title{Implementation of SVM plus}

\description{
Implementation of SVM plus classification problems. See (Vladimir et. al, 2009, <https://doi.org/10.1016/j.neunet.2009.06.042>) for theoretical details and see (Li et. al, 2016, <https://github.com/okbalefthanded/svmplus_matlab>) for implementation in MATLAB.
}


\details{
LUPI is a novel machine learning paradigm. It offers faster convergence of the learning process by exploiting the privileged information. In other words, ``fewer training examples
are needed to achieve similar predictive performance" or ``the
same number of examples can provide a better predictive performance".	In SVM+, privileged  information is used to estimate a linear model of the slack variables, namely \eqn{$\xi_i = (w^*)^T x_i^* + b^*$}.

The SVM+ objective function is defined as:
\eqn{
 \min_{w,b} && \frac{1}{2} w^T w \color{red}{+\frac{\gamma}{2} (w^*)^T (w^*) + C \sum_{i=1}^l [(w^*)^T x_i^* + b^*]} }
 \eqn{
s.t. && y_i (w^T x_i + b) \geq 1- \color{red}{[(w^*)^T x_i^* + b^*], \; (w^*)^T x_i^* + b^* \geq 0}, \;  \color{black}\forall i
}
}
