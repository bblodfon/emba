% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_link_operator_diff_mat_based_on_mcc_clustering}
\alias{get_avg_link_operator_diff_mat_based_on_mcc_clustering}
\title{Get average link operator difference matrix based on MCC clustering}
\usage{
get_avg_link_operator_diff_mat_based_on_mcc_clustering(
  models.mcc,
  models.link.operator,
  num.of.mcc.classes,
  penalty = 0
)
}
\arguments{
\item{models.mcc}{a numeric vector of Matthews Correlation Coefficient (MCC)
scores, one for each model. The \emph{names} attribute holds the models' names.
Can be the result of using the function \code{\link{calculate_models_mcc}}.}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes.
The row names specify the models' names (same order as in the
\code{models.mcc} parameter) whereas the column names specify the
network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
(\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
both activating and inhibiting regulators (no link operator).}

\item{num.of.mcc.classes}{numeric. A positive integer larger than 2 that
signifies the number of mcc classes (groups) that we should split the models
MCC values.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a matrix whose rows are \strong{vectors of average node link operator
differences} between two groups of models where
the classification was based on the models' MCC values.
Rows represent the different classification group matchings, e.g. (1,2) means the
models that belonged to the 1st group of MCC values vs the models that
belonged to the 2nd group. The columns represent the network's node names.
Values are in the [-1,1] interval.
}
\description{
This function uses the \code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}}
function with the parameter \code{models.link.operator} as input in the place of
\code{models.stable.state}, since the two matrices representing the two inputs
have the same data format (rows represent models, columns represent nodes,
and each value is a number in the [0,1] interval).
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node's boolean equation has the \strong{AND NOT} link operator in the
'good' models compared to the 'bad' ones while a value closer to 1 means that
the node's boolean equation has mostly the \strong{OR NOT} link operator
in the 'good' models. A value closer to 0 indicates that the link operator in
the node's boolean equation is \strong{not so much different} between the
'good' and 'bad' models and so it won't not be a node of interest when
searching for indicators of better performance (higher average MCC value)
in the parameterization of the good models (the boolean equations). A value
exactly equal to 0 can also mean that this node didn't not have a link operator
in its boolean equation, again making it a non-important indicator of difference
in model performance.
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
