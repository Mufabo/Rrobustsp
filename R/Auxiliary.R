#' normalized median absolute deviation
#'
#' madn computes the normalized median absolute deviation estimate of
#' scale, i.e.,
#'
#'mu_hat = arg min_mu SUM_i rho_TUK(y_i - mu)
#'
#'
#'@param y : data vector of size N x 1
#'
#'
#'@return sig : normalized median absolute deviations scale estimate
#'
#'@examples
#'madn(rnorm(5))
#'
#'@importFrom stats median
#'
#'@export
madn <- function(y)
    {
    # Compute the auxiliary scale estimate as
    if (is.complex(y)) const <- 1.20112 else const <- 1.4815
    return (const*median(abs(y-median(y))))
}



#'Huber's score function
#'
#'@param x is N x 1 data vector x which can be complex
#'or real
#'
#'@param threshold constant c
#'@export
psihub <- function(x,c){
    if(is.complex(x)){
        signum_x <- x / abs(x)
        return (x*(abs(x)<=c) + c*signum_x * (abs(x)>c))
    }
    else{return (x*(abs(x)<=c) + c*sign(x)*(abs(x)>c))}
}

#'@export
psituk <- function(x, c){
    return (x*((1-(abs(x)/c)^2)^2 )*(abs(x)<=c))
}

#'@export
scaledata <- function(x){
    nom <- sweep(x, 2, apply(x, 2, min)) # subtract col_mins from each respective column
    denom <- apply(x, 2, max) - apply(x, 2, min)
    return(3 * sweep(nom, 2, denom, FUN = '/'))
}

#'@export
rhohub <- function(x, c){
    return ((abs(x)^2)*(abs(x)<=c) + (2*c*abs(x)-c^2)*(abs(x)>c))
}

#'@export
rhotuk <- function(x, c){
    return ((c^2/3)*((1-(1-(abs(x)/c)^2)^3)*(abs(x)<=c) + (abs(x)>c) ))
}

#'@export
soft_thresh <- function(x, t){
    s <- abs(x) - t
    s <- (s + abs(s)) / 2
    return (if (is.complex(x)) x*s/abs(x) else sign(x)*s)
}

#'@export
whub <- function(absx, c){
    wx = 1*(Re(absx)<=c) + (c*(1/absx))*(Re(absx)>c)
    wx[absx==0] = 1
    return (wx)
}

#'@export
wtuk <- function(absx, c){
    return (((1-(absx/c)^2)^2)*(Re(absx)<=c))
}
