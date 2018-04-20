\docType{package}
\name{svmplus}
\alias{svmplus}
\title{Implementation of SVM Plus}

\description{
Implementation of SVM plus for classification problems.
}


\details{
The classical machine learning paradigm assumes, training examples in the form of iid pair:
\eqn{$(x_1, y_1), ..., (x_l, y_l), \hspace{1em} x_i \in X, \hspace{1em} y_i \in \{-1, +1\}.$}

Training examples are represented as features \eqn{x_i} and the same feature space is required for predicting future observations. However, this approach does not make use of other useful data that is only available at training time; such data is referred to as Privileged Information (PI).

Learning Under Privileged Information (LUPI) is a novel machine learning paradigm. It offers faster convergence of the learning process by exploiting the privileged information. In other words, ``fewer training examples are needed to achieve similar predictive performance" or ``the
same number of examples can provide a better predictive performance".	In LUPI paradigm, training examples come in the form of iid triplets

\eqn{$(x_1,x_1^*, y_1), ..., (x_l,x_l^*, y_l), \hspace{1em} x_i \in X, \hspace{1em} x^*_i \in X^*, \hspace{1em} y_i \in \{-1, +1\} $}

where \eqn{x^*} denotes PI. SVM+ is one realization of LUPI paradigm.
In SVM+, privileged  information is used to estimate a linear model of the slack variables, namely


\eqn{$\xi_i = (w^*)^T z_i^* + b^*,$}

 where \eqn{ z_i = \phi(x_i)} represents the kernel mapping.


The SVM+ objective function is defined as:
\eqn{$ \min_{w,b} \left\lbrace \frac{1}{2} w^T w +\frac{\gamma}{2} (w^*)^T (w^*)
+ C \sum_{i=1}^l [(w^*)^T z_i^* + b^*] \right\rbrace $}


 \eqn{$ s. t. \quad y_i (w^T z_i + b) \geq 1- [(w^*)^T z_i^* + b^*], $}


 \eqn{$ (w^*)^T z_i^* + b^* \geq 0,   \forall i$}



The dual SVM+ problem is defined as follow.

\eqn{$ \max_{w,b} \left\lbrace \sum_{i=1}^l \alpha_i - \frac{1}{2} \sum_{i,j=1}^l \alpha_i \alpha_j
y_i y_j K(x_i, x_j) -  \frac{1}{2\gamma} \sum_{i,j=1}^l (\alpha_i+\beta_i - C) (\alpha_j+\beta_j - C) K^*(x_i^*, x_j^*)
\right\rbrace $}


 \eqn{$s. t. \quad \sum_{i=1}^l \alpha_i y_i = 0, \quad \sum_{i=1}^l (\alpha_i+\beta_i - C) = 0, $}


\eqn{$ \alpha_i \geq 0, \quad \beta_i \geq 0 $}



This package offeres a Quadratic Programming (QP) based convex optimization solution for the dual SVM+ problem. In future, LIBSVM and LibLinear based faster implementaions are planned to be supported. We refer to [1] for theoretical details of LUPI and SVM+, and we refer to [2] for implementation details of SVM+ in MATLAB.


}

\references{
[1] Vladimir et. al, Neural Networks, 2009, 22, pp 544–557.
\url{https://doi.org/10.1016/j.neunet.2009.06.042}

[2] Li et. al, 2016.
\url{https://github.com/okbalefthanded/svmplus_matlab}

[3] Bendtsen, C., et al., Ann Math Artif Intell, 2017, 81, pp 155–166.
\url{https://doi.org/10.1007/s10472-017-9541-2}
}
