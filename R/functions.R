## A is the abundance matrix and T is the trait matrix and dset is the data set

# Community weighted means for each modality within each bromeliad -> we need to standardize for every bromeliad 
# Fuzzit produces for every bromeliad the traits weighted by the abundance of every species in a plant. 


Fuzzit4a<-function(A,T)
{
  region.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			# bromeliad X traits
  n.plant<-dim(A)[2]
#  A<-sqrt(A)
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums(sqrt(A[,i]/sum(A[,i]))*T,na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(region.id,a.fuz)
  return(a.fuz3)
}


permut.row.matrix <- function (data, strata = NULL, seqpermutation = NULL) 
{
  N <- dim(data)[1]
  if (!is.null(strata) & N != length(strata)) {
    stop("\n strata must be the length of number of row in the data\n")
  }
  if (!is.null(seqpermutation) & N != length(seqpermutation)) {
    stop("\n seqpermutation must be the length of number of row in the data\n")
  }
  if (is.null(seqpermutation)) {
    samp <- permut.vector(N, strata = strata, nset = 1)
  }
  else {
    samp <- seqpermutation
  }
  permut.matrix <- data[samp, , drop = FALSE]
  res <- list(permut.matrix = permut.matrix, samp = samp)
  return(res)
}