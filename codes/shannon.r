# Return the Shannon entropy & Shannon variation.
# A : array of read-numbers.
shannon_entropy <- function(A) {
    T <- 0
    S <- 0
	for (a in A) {
        if (a > 0) {
            T <- T + a
            S <- S + log(a)*a
        }
	}
    return(-S/T + log(T))
}

shannon_variation <- function(A) {
    return(exp(shannon_entropy(A)))
}

V <- c(1,1,1,5)
print(shannon_variation(V))
