\name{plot.ivc}
\alias{plot.ivc}
\alias{ivc_plot}
\title{Plot the IVC diagnostic}
\description{
Draws a forest plot with all three quantities on a common scale, each with
its own confidence interval: the conventional IV estimate, the
compass-corrected estimate, and their contrast Delta. When the interval for
Delta excludes zero, that row is shaded and colored to flag evidence against
exogeneity. A base-graphics plot is used so the package has no plotting
dependencies.
}
\usage{
\method{plot}{ivc}(x, ...)

ivc_plot(x, ...)
}
\arguments{
  \item{x}{An object of class \code{"ivc"} from \code{\link{compare_iv_cv}}.}
  \item{...}{Further arguments passed to the underlying plotting calls.}
}
\value{
\code{x}, invisibly. Called for the side effect of drawing a plot.
}
\examples{
d <- ivc_simulate(n = 1000, gamma = 1, kappa = 0.3, seed = 1)
fit <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")
plot(fit)
}
