## A is the abundance matrix and T is the trait matrix and dset is the data set

# Community weighted means for each modality within each bromeliad -> we need to standardize for every bromeliad 
# Fuzzit produces for every bromeliad the traits weighted by the abundance of every species in a plant. 

# THINGS TO CHANGE
# W is proportional before multiplies B
# W needs to be divided by the total number of individuals in each bromeliad
# Proportional representation of traits

Fuzzit<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])	#bromeliad x  mean fuzzy coded trait matrix
  n.plant<-dim(A)[2]
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums(A[,i]*T,na.rm=TRUE)/n.plant		# average [A]abundance weighted traits for a plant
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz2<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz2)[2]<-"dset"
  return(a.fuz2)
}


##Weight by the proportional abundance

Fuzzit2<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			
  n.plant<-dim(A)[2]
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums(A[,i]/sum(A[,i])*T,na.rm=TRUE)		# proportional [A]bundance weighted [T]raits for a plant
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz2<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz2)[2]<-"dset"
  return(a.fuz2)
}

##Hellinger transformation of abundance

Fuzzit3<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			
  n.plant<-dim(A)[2]
  A<-sqrt(A)
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums((A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz3)[2]<-"dset"
  return(a.fuz3)
}

Fuzzit3a<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			
  n.plant<-dim(A)[2]
  #  A<-sqrt(A)
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums(sqrt(A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz3)[2]<-"dset"
  return(a.fuz3)
}

Fuzzit3b<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			
  n.plant<-dim(A)[2]
  A<-sqrt(A)
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums((A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz3)[2]<-"dset"
  return(a.fuz3)
}

Fuzzit3c<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			
  n.plant<-dim(A)[2]
  #  A<-sqrt(A)
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums(log(A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz3)[2]<-"dset"
  return(a.fuz3)
}

Fuzzit3d<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			
  n.plant<-dim(A)[2]
  A<-log(A+1)
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums((A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz3)[2]<-"dset"
  return(a.fuz3)
}

#Kurt's version
Fuzzit4<-function(A,T)
{
  region.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			# bromeliad X traits
  n.plant<-dim(A)[2]
  A<-sqrt(A)
  for(i in 1:dim(A)[2]){
    #    a.fuz[i,]<-colSums(sqrt(A[,i]/sum(A[,i]))*T,na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
    a.fuz[i,]<-colSums((A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(region.id,a.fuz)
  return(a.fuz3)
}

Fuzzit4a<-function(A,T)
{
  region.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			# bromeliad X traits
  n.plant<-dim(A)[2]
#  A<-sqrt(A)
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums(sqrt(A[,i]/sum(A[,i]))*T,na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
#    a.fuz[i,]<-colSums((A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz3<-data.frame(region.id,a.fuz)
  return(a.fuz3)
}

Fuzzit4b<-function(A,T)
{
  region.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])			# bromeliad X traits
  n.plant<-dim(A)[2]
  A<-sqrt(A)
  for(i in 1:dim(A)[2]){
    #    a.fuz[i,]<-colSums(sqrt(A[,i]/sum(A[,i]))*T,na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
    a.fuz[i,]<-colSums((A[,i]/sum(A[,i])*T),na.rm=TRUE)		# square-root transforme proportional [A]bundance weighted [T]raits for a plant  
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