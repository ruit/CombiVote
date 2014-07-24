#Tian R. <tianremi@gmail.com>
#July 23, 2014
#each P(AB|C) is a weak classifier
#sample(c("posi", "nega"), 1, prob=c(P(AB|"posi"), P(AB|"nega"))) 
#might be useful in cases where features show interactions
#sample partial feature pairs
#divide training data into 3 equal parts

#weights: random 1 or 0
#suppose 100 feature pairs
#w<-sample(c(0,1), 100, replace=T, prob=c(0.2, 0.8))
#mutate one position (0->1, or 1->0) each time
#if perf improves, keep this mutation; otherwise, keep as it was
#move onto to next positions, in total 100 positions

#Invention of a new tool is not easy!!!

#input is df of dicretized data

#@July 24, 2014
# sample part of the features each time to avoid overfitting
# K-means clustering to group/ discretize data based on 2 or more features
# boostrap to get better estimation of conditional probabilities

Fill.NA <- function(mat, method=c("median", "mean")){

    method <- match.arg(method)
    
    if(method == "median"){
        output <-
            apply(mat, 2, function(x) x[is.na(x)] <- median(x, na.rm=TRUE))
    }
:w
    else if(method == "mean"){
        output <-
            apply(mat, 2, function(x) x[is.na(x)] <- mean(x, na.rm=TRUE))
    }
    return (output)
}


combiVote<-function(data, k=ncol(data), epsilon=0.001, ...){
	#check input
	if (is.matrix(data)) {
		data<-as.data.frame(data)
		cat ("Ok!\n")
	} else  if (is.data.frame(data)){
		cat ("good!\n")
	} else {
		stop ("Input must be either array or dataframe!")
				}	
	#partition data into 3 parts, train1, train2, test
	labels<-names(table(data[,k]))
	train1_index<-c()
	for (i in 1:length(labels)){
		num<-nrow(data[data[,k]==labels[i],])
		index<-which(data[,k]==labels[i])
		train1_index<-c(train1_index,
			sample(index, round(num/3, 0)))	
		
		}
	train1Data<-data[train1_index,]
	rest<-data[-train1_index,]

	train2_index<-sample(1:nrow(rest), round(0.5*nrow(rest), 0))
	train2Data<-rest[train2_index,]

	testData<-rest[-train2_index,]
	###################check above!##################
	
	#train for the 1st round

	modelOne<-naiveBayes(train1Data[,-k], train1Data[,k])	
	ori.model<-modelOne

	biaoxian<-function (alpha, modelOne=modelOne){
		for (i in 1:length(modelOne$tables)){
                        modelNew$tables[[i]]<-alpha[i]*modelOne$tables[[i]]###for discrete data only!!!!
                        }
		modelOne<-modelNew
		ptable<-table(predict(modelOne, train2Data[,-k], type="class"),train2Data[,k])
		p1<-ptable[1, 1]/sum(ptable[,1])
		p2<-ptable[2,2]/sum(ptable[,2]) 
		perf1<-0.5*(p1+p2)/((p1-p2)^2+1)	
		return (perf1)
		}
	
	
	#train for the 2nd round 
	# assign a stochastic weight to each feature sequentially using uniform distribution
		result <- optimize(biaoxian, c(0,1))
		
		#updated modelOne
		#orginal modeloriginal
		#w<-result$alpha
		#opt.model<-ori.model
		#for (i in 1:length(modelOne$tables)){
                 #       opt.model$tables[[i]]<-w[i]*ori.model$tables[[i]]###for discrete data only!!!!
                 #       }
		return (result)
	}
