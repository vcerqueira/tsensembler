% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build-base-ensemble.r
\name{learning_base_models}
\alias{learning_base_models}
\title{Training the base models of an ensemble}
\usage{
learning_base_models(train, form, specs, num_cores)
}
\arguments{
\item{train}{training set to build the predictive models;}

\item{form}{formula;}

\item{specs}{object of class \code{\link{model_specs-class}}}

\item{num_cores}{A numeric value to specify the number of cores used to
train base and meta models. num_cores = 1
leads to sequential training of models. num_cores > 1
splits the training of the base models across num_cores cores.}
}
\value{
A series of predictive models (\code{base_model}), and
the weights of the models computed in the training
data (\code{preweights}).
}
\description{
This function uses \emph{train} to build a set
of predictive models, according to \emph{specs}
}
\examples{
data("water_consumption")
dataset <- embed_timeseries(water_consumption, 5)
specs <- model_specs(c("bm_ppr","bm_svr"), NULL)
M <- build_base_ensemble(target ~., dataset, specs, 1)

}
\seealso{
\code{\link{build_base_ensemble}}.
}
