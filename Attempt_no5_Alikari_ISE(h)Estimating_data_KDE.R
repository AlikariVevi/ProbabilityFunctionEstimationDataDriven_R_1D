
#############################################################################
#############################################################################
#######################        DATA SAMPlE        ###########################
#############################################################################
#############################################################################

mainLabels<-c("Significant Wave height, Iraklio","Wave Energy Period, Iraklio","Mean Direction of Wave, Iraklio")
xLabels<-c("Significant Wave height (m)", "Wave Energy Period (sec)", "Mean Direction of Wave")

#getwd()

mydata<-read.table(file="data.txt")

####################################
########    Significant    #########
########    Wave height    #########
####################################

#u<-mydata[,5]
#dataprint<-print("Significant Wave height, Iraklio")
#mainLabel<-mainLabels[1]
#xLabel<-xLabels[1]

#file_name="SUMMARY_Significant_Wave_height_Iraklio"



####################################
########    Wave Energy    #########
########       Period      #########
####################################

#u<-mydata[,6]
#dataprint<-print("Wave Energy Period, Iraklio")
#mainLabel<-mainLabels[2]
#xLabel<-xLabels[2]

#file_name="SUMMARY_Wave_Energy_Period_Iraklio"


####################################
#######    Mean Direction   ########
#######      of Wave        ########
####################################

u<-mydata[,7]
dataprint<-print("Mean Direction of Wave, Iraklio")
mainLabel<-mainLabels[3]
xLabel<-xLabels[3]

file_name="SUMMARY_Mean_Direction_of_Wave_Iraklio"


##############################################################################
##############################################################################
##############################################################################
##############################################################################

file_variable=file(file_name, "w")
nprint<-sprintf("number of data n=%i",length(u))
write(nprint,file=file_name,append = TRUE)

#############################################################################
#############################################################################
#######################      KERNEL FUNCTIONS     ###########################
#############################################################################
#############################################################################


####################################
########       Gauss       #########
########  Kernel Function  #########
####################################


K<-function(u){
return((sqrt(2*pi))^(-1)*exp(-((u)^2)/(2)))
}

Kernelprint<-print("Kernel function : Gauss")
write("\nKernel function : Gauss\n",file_name,append=TRUE)

####################################
########      Uniform      #########
########  Kernel Function  #########
####################################

#K<-function(u){
#return(1/2)
#}

#Kernelprint<-print("Kernel function : Uniform")
#write("\nKernel function : Uniform\n",file_name,append=TRUE)

####################################
########   Epanechnikov    #########
########  Kernel Function  #########
####################################

#K<-function(u){
#return((3/4)*(1-u^2))
#}

#Kernelprint<-print("Kernel function : Epanechnikov")
#write("\nKernel function : Epanechnikov\n",file_name,append=TRUE)


####################################
########     Silverman     #########
########  Kernel Function  #########
####################################

#K<-function(u){
#return((1/2)*exp(-abs(u)/sqrt(2))*sin((abs(u)/sqrt(2))+(pi/4)))
#}

#Kernelprint<-print("Kernel function : Silverman")
#write("\nKernel function : Silverman\n",file_name,append=TRUE)


##############################################################################
##############################################################################
##############################################################################
##############################################################################

####################################
########      Scaled       #########
########  Kernel Function  #########
####################################

Kh<-function(h,u){
return((1/h)*K(u/h))
}

####################################
########   Kernel Density   ######## 
######## Estimator Function ########
####################################

f<-function(x,u,h){
n=length(u)
return((1/n)*sum(Kh(h,x-u)))
}

####################################
###### Histograms parameters #######
####################################

####################################
############# breaks ###############
####################################

  ## the number of breaks 
  ##       effect the
  ## estimated density fanction 
  ## by effecting the h parameter

###########################################
############ Sturges' Rule ################
###########################################

  ## Sturges' rule considers the binomial 
  ## distribution as a model of an 
  ## optimally costructed histogram
  ## this rule is an estimation of 
  ## the number of breaks

	## the smallest integers not less than the
	## the computed value
#br=ceiling(1+log2(length(u)))

  ## build in sturges' rule 
#br=nclass.Sturges(u)

#br_methode_print<-sprintf("histogram's optimally construction rule : Sturges\n")
#hist(u, prob=TRUE,breaks=br,main="Struges Rule",xlab=xLabel)

###########################################
##### Normal Density Reference Rule #######
#####             Scott             #######
###########################################

  ## normal density reference rule 
  ## uses the samples standar deviation
  ## to optain an optimall histogram
  ## this rule is an estimation of 
  ## the width of the breaks

#br_width=3.5*sd(u)*length(u)^(-1/3)
	## the smallest integers not less than the
	## the computed value
#br=ceiling((max(u)-min(u))/br_width)

  ## build in scott's rule 
#br=nclass.scott(u)

#br_methode_print<-sprintf("histogram's optimally construction rule : Scott\n")
#hist(u, prob=TRUE,breaks=br,main="Scott's Rule",xlab=xLabel)

###########################################
#####    Freedman-Diaconis Rule     #######
###########################################

  ## Freedman-Diaconis rule 
  ## 
  ## 
  ## this rule is an estimation of 
  ## the width of the breaks



  ## build in Freedman-Diaconis's rule 
br=nclass.FD(u)

br_methode_print<-sprintf("histogram's optimally construction rule : Freedman-Diaconis\n")
hist(u, prob=TRUE,breaks=br,main="Freedman Diaconis Rule",xlab=xLabel)

################################################
########### print and write result #############
################################################

write(br_methode_print,file=file_name,append = TRUE)
brprint<-sprintf("histogram's breaks=%f",br)
write(brprint,file=file_name,append = TRUE)

#############################################################################
#############################################################################
#######################          BANDWIDTH        ###########################
#############################################################################
#############################################################################

####################################
##### Overestimated Value of #######
##### of smoothing parameter #######
####################################

h=10

####################################
####################################
######   Finding optimall   ########
######  smoothing parameter ########
####################################
####################################

####################################
#### BEGIN h ESTIMATING FUNCTION ###
####    summed squared error     ###
####################################


SqErrorEstimating<-function(h,f,u,br){
#########  histogram  ##############
	hi<-hist(u,breaks=br)
#######  histogram's mids  #########
	midV=hi$mids
#######  density histogram  ########
	dens=hi$density
######  Applying estimator  ########
###### on mids of histogram ########
	datafmid<-sapply(midV,f,u,h)
######    squared error    #########
	ISEh<-0
	for (i in 1:length(midV)){
		ISEh=ISEh+(datafmid[i]-dens[i])^2
		}
	return(ISEh)
}

####################################
#### 	   estimating minimum     ####
####    summed squared error    ####
####################################
#e<-10^(-3.1)

###################################
minimumError<-function(h1,h2,f,u,br){
Value1<-SqErrorEstimating(h1,f,u,br)
Value2<-SqErrorEstimating(h2,f,u,br)
if(Value1<Value2){
return (Value1)
}else {
return (Value2)
}
}

htominimumError<-function(h1,h2,f,u,br){
Value1<-SqErrorEstimating(h1,f,u,br)
Value2<-SqErrorEstimating(h2,f,u,br)
if(Value1<Value2){
return(h1)
}else {
return(h2)
}
}

###################################
#### END h ESTIMATING FUNCTION ####
####  integrated squared error ####
###################################


### initializing checkinh loop ###

hBoundary<-c(h,h)

#### checking loop fo h #####
nn=50
for (i in 1:nn){
ISEh<-minimumError(hBoundary[1],hBoundary[2],f,u,br)
h<-htominimumError(hBoundary[1],hBoundary[2],f,u,br)
hBoundary[1]<-h
hBoundary[2]<-h-h/2
}



######### printable strings ###########

hprint<-sprintf("Estimated bandwidth h=%f",h)
write("\n",file=file_name,append=TRUE)
write(hprint,file=file_name,append=TRUE)


hestimationPrint1<- sprintf("the smoothing parameter has been estimated")
hestimationPrint2<- sprintf("by assuming that histogram's densities are")
hestimationPrint3<- sprintf("the true density values for the middle of each")
hestimationPrint4<- sprintf("bar.")
hestimationPrint5<- sprintf("Thus a squared error sum over these points")
hestimationPrint6<- sprintf("is been used to accomplice the desirable.")
hestimationPrint7<- sprintf("Squared estimated error=%f",ISEh)


##############################################################################
##############################################################################
##############################################################################
##############################################################################



####################################
####################################
####    Applying estimator     #####
####          on data          #####
####       with the new        #####
####    smoothing parameter    #####
#### on data (training inputs) #####
####################################
####################################

dataf<-sapply(u,f,u,h)


#############################################################################
#############################################################################
#######################           GRAPH           ###########################
#############################################################################
#############################################################################

####################################
#########     Histogram     ########
####################################

## for changing y-axes add ylim=c(min,max)
hist(u, prob=TRUE,breaks=br,main=mainLabel,xlab=xLabel)

####################################
#########        KDE        ########
####################################

######### training data ########
#points(u,dataf,pch=20, cex=0.02, col="red")
#legend(1.5,0.7,"estimated density on data",lty=3,col="black",bty="n")

rank<-order(u)
lines(u[rank],dataf[rank], type="l", col="red", lwd=2)
legend(legends_coord <- locator(1) ,"Vevi's Algorithm",lty=1,col="red",bty="n")


#############################################################################
#############################################################################
#######################           SUMMARY           #########################
#############################################################################
#############################################################################

nprint
dataprint
Kernelprint
br_methode_print
brprint
hprint
hestimationPrint1
hestimationPrint2
hestimationPrint3
hestimationPrint4
hestimationPrint5
hestimationPrint6




close(file_variable)




#############################################################################
#############################################################################
#################     R kernel estimation packages    #######################
#############################################################################
#############################################################################

####################################
####    KernelSmooth package    ####
#### Wand, M.P. and Jones, M.C. ####
####################################

#install.packages("KernSmooth")
library("KernSmooth")

EstimationWand<-bkde(u)

lines(EstimationWand,main=mainLabel,xlab=xLabel,col="blue",type="l", lwd=2)
legend(legends_coord <- locator(1) ,"R KernSmooth-method",lty=1,col="blue",bty="n")

####################################
#########   ks package   ###########
#########   Tarn Duong   ###########
####################################

#install.packages('ks')
library('ks')

EstimationDuong<-kde(u)


plot(EstimationDuong,main=mainLabel,xlab=xLabel,col="green",lty=3,lwd=4,add=TRUE)
legend(legends_coord <- locator(1) ,"R ks-method",lty=3,col="green",bty="n")

