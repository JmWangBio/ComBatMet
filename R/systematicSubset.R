systematicSubset <- function(n,order.by)
#	Copied from edgeR
{
	ntotal <- length(order.by)
	sampling.ratio <- floor(ntotal/n)
	if(sampling.ratio <= 1) return(1:ntotal)
	i1 <- floor(sampling.ratio/2)+1
	i <- seq.int(from=i1,to=ntotal,by=sampling.ratio)
	o <- order(order.by)
	o[i]
}
