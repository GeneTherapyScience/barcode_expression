# Return the smallest read-number to keep.
# A : array of read-numbers.
wsn_thres <- function(A, ratio=10^(-5)) {
	A <- sort(A)
	S <- sum(A)
	prev <- 0
	for (a in A) {
		if (a == prev) {
			S <- S - a
		} else if (a < S*ratio) {
			S <- S - a
			prev <- a
		} else {
			return(a)
		}
	}
	return(a)
}

V <- c(1,4,1,5,25,20,1,1,1,3,1,2,2,2,1,1)
# V = c(1,1,1,1)
r <- 0.1
print(wsn_thres(V,r))
