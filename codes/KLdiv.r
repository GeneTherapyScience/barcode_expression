# Kullback-Leibler divergence
KLdiv <- function(P, Q) {
    S <- 0
    for (i in 1:length(Q)) {
        if (P[i] > 0) {
            S <- S + P[i]*(log(P[i])-log(Q[i]))
        }
    }
    return(S)
}

# Jensen-Shannon divergence
JSdiv <- function(A, B) {
    M <- (A+B)/2
    return((KLdiv(A,M) + KLdiv(B,M))/2)
}

A <- c(1/2,1/4,1/4,0)
B <- c(0,1/4,1/2,1/4)
print(JSdiv(A,B))
