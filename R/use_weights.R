#' Get Importance Resampling Indices From Weights
#'
#' Takes a vector of weights, or data frame or list containing sets of weights,
#' and resamples indices for use in later computation.
#'
#' @param x A vector of weights, a list of weight vectors, or a data frame of
#'   type \code{adjustr_weighted} containing a \code{.weights} list-column
#'   of weights.
#' @param frac A real number giving the fraction of draws to resample; the
#'   default, 1, resamples all draws. Smaller values should be used when
#'   \code{replace=FALSE}.
#' @param replace Whether sampling should be with replacement. When weights
#'   are extreme it may make sense to use \code{replace=FALSE}, but accuracy
#'   is not guaranteed in these cases.
#'
#' @return A vector, list, or data frame, depending of the type of \code{x},
#' containing the sampled indices. If any weights are \code{NA}, the indices
#' will also be \code{NA}.
#'
#' @export
get_resampling_idxs = function(x, frac=1, replace=T) {
    if (frac < 0) stop("`frac` parameter must be nonnegative")
    get_idxs = function(w) {
        if (all(is.na(w))) return(NA_integer_)
        sample.int(length(w), size=round(frac*length(w)), replace=replace, prob=w)
    }

    if (is(x, "list")) {
        map(x, get_idxs)
    } else if (is(x, "adjustr_weighted")) {
        x$.idxs = map(x$.weights, get_idxs)
        x
    } else {
        get_idxs(x)
    }
}
