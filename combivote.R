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

#@July 25, 2014, keep strong, keep fucking Japs!

Fill.NA <- function(mat, method=c("median", "mean")){
#what about knn?
    method <- match.arg(method)
    
    if(method == "median"){
        output <-
            apply(mat, 2, function(x) x[is.na(x)] <- median(x, na.rm=TRUE))
    }

    else if(method == "mean"){
        output <-
            apply(mat, 2, function(x) x[is.na(x)] <- mean(x, na.rm=TRUE))
    }
    return (output)
}



Fscore<-function (x1, x2, na.rm=T, epsilon=0.00001){

    #x1, x2 are two numeric vectors
    if (!is.numeric(x1) | !is.numeric(x2) ) 
    stop ("Input of Fscore should be 2 numeric vectors!") 
    
    if (na.rm){
        p<-x1[!is.na(x1)]
	n<-x2[!is.na(x2)]
        }
    
    comb<-c(p, n)
    
    #scale the data 
    posi<-(p-mean(comb))/sd(comb)
    nega<-(n-mean(comb))/sd(comb)    
 
    total_mean<-mean(c(posi, nega))
    posi_mean<-mean(posi)
    nega_mean<-mean(nega)
    
    numerator<-(posi_mean-total_mean)^2+(nega_mean-total_mean)^2
    denominator<-var(posi)+var(nega)

    if (denominator == 0){
        denominator<-denominator+epsilon}
    
    fscore<-round(sqrt(numerator/denominator), 3)
    
    return (fscore)
	
    }
 


RankFeaByImportance<-function (aml2="data.frame", k=ncol(aml2), top=20, ...){
    
    #Dependant on Fscore function, k is the col for class (labels)
    #aml2 needs to be scaled
    if (!is.data.frame(aml2)) stop ("Input should be a data frame!")
    
    colnames(aml2)[k]<-"class"
    tb<-table(aml2$class)
    
    #consider binnary classification problems
    zheng<-attr(tb,"names")[1]
    fu<-attr(tb,"names")[2]
	
    fscore<-c()
    for (i in (1:ncol(aml2))[-k]){    
        fscore<-c(fscore,Fscore(aml2[aml2$class==zheng,i], 
                                aml2[aml2$class==fu, i])
                 )
	}
    names(fscore)<-colnames(aml2)[-k]
    
    ##barplot(head(sort(fscore, decreasing=T), 50), las=2, cex.names=0.75, main="Top 50 most discriminating protein features")
    fea<-data.frame(fscore=sort(fscore, decreasing=T))
    
    index<-c()
    for (i in 1:top){
        index<-c(index, which(colnames(aml2)==rownames(fea)[i]))
        }
    
    cat ("The F scores of the most important features are:\n")	
    print (fea[1:top,])

    return (aml2[,c(k,index)])
    }



Disc.quantile<-function(vect){

    #
    if (!is.numeric(vect)) stop ("Input of function Disc.quantile must be a numeric vector!")
    
    new<-rep("NULL", length(vect))

    new[vect<quantile(vect,probs=0.25)]<-"L1"
    new[vect>=quantile(vect,probs=0.25)]<-"L2"
    new[vect > quantile(vect,probs=0.75)]<-"L3"
    return(new)
    }



Get_Kcut<-function(HMdnt.vector, k=4)
{	
    breakpoints<-c()
    brk.mat<-c()
    for (N in 1:100){
	    fit <- kmeans(HMdnt.vector, k)
		aggregate(HMdnt.vector,by=list(fit$cluster),FUN=mean)
		cell.kmeans <- data.frame(HMdnt.vector, fit$cluster)
	    colnames(cell.kmeans)<-c("histM", "Cst")
#table(cell.kmeans$cluster)
#Figure out the order of "1 2 3 4"as assigned by kmeans
        for (i in 1:k){
	        cluster<-cell.kmeans[cell.kmeans$Cst==i, ]
			breakpoints<-c(breakpoints, min(cluster[,1]), max(cluster[,1]))
		}
		brk.s<-sort(breakpoints)
	    breakpoints<-c()
		brk.mat<-rbind(brk.mat,brk.s)
	}
	
    cut<-c()
	for (i in 1:(2*k)){
#   	k=4, breakpoints=8
	    cut<-c(cut, median(brk.mat[,i]))
	}
	return (cut)
}

############################
Get_level<-function(x, vector){#vector is the cutting points
	x.lev<-c()
	for (i in 1:length(x)){
		if (x[i]>=vector[1] && x[i]<=vector[2]) {
			x.lev[i]<-"Lev0"
		} else {
			if (x[i]>=vector[3] && x[i]<=vector[4]){
				x.lev[i]<-"Lev1" 
			}  else {
				if(x[i]>=vector[5] && x[i]<=vector[6]){
					x.lev[i]<-"Lev2" 
				} else {
					if(x[i]>=vector[7] && x[i]<=vector[8]){x.lev[i]<-"Lev3" }
				}
			}
		}    
		
	}
	return (x.lev)
	
}



combiVote<-function(data, k=ncol(data), epsilon=0.001, ...){

	#check input
	if (is.matrix(data)) {
		data<-as.data.frame(data)
		#cat ("Ok!\n")
	} else  if (is.data.frame(data)){
		#cat ("good!\n")
	} else {
		stop ("Input must be either matrix or dataframe!")
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
