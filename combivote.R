#Tian R. <tianremi@gmail.com>
#July 23, 2014

#require glmnet
#familiy="Binomial" or "Multinomial"
#alpha=1, means L1 lasso penalty

#each P(AB|C) is a weak classifier
#sample(c("posi", "nega"), 1, prob=c(P(AB|"posi"), P(AB|"nega"))) 
#might be useful in cases where features show interactions

#weights: random 1 or 0

#suppose 100 feature pairs
#w<-sample(c(0,1), 100, replace=T, prob=c(0.2, 0.8))
#mutate one position (0->1, or 1->0) each time
#if perf improves, keep this mutation; otherwise, keep as it was
#move onto to next positions, in total 100 positions


#@July 24, 2014
#@July 25, 2014, keep strong, keep fucking Japs!
#@July 29, 2014, keep calm, ensure safety.
#@July 30, 2014, make every good hour, each day, week!


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



Disc.kmeans<-function(mat, k=9, times=100){	
    
    #mat should be scaled
    if (!is.matrix(mat)) stop ("The input of Disc.kmeans must be matrix!")
    mat<-scale(mat)

    output<-c()
    for (N in 1:times){
        fit<-kmeans(mat, k)
	
	#zhongxin has 4 rows, number of columns is ncol(mat)
	zhongxin<-fit$centers
	#print (zhongxin)
	
	zhongxin<-zhongxin[order(zhongxin[,1],zhongxin[,2]),]#depends on the ncol (mat)
       	#cat ("Sorted centers are:\n") 
        #print (zhongxin)
        
        #function to assign a point to its nearest centroid
	group<-function(vect){

            #dist() is a class not a mtrix
            juli<-as.matrix(dist(rbind(vect, zhongxin), method="euclidean"))[-1,1]
            g<-which(order(juli)==1)
	    return (paste("L", as.character(g), sep=""))
            }
        
        output<-cbind(output,  apply(mat, 1, group))

        }  	

    major<-function (x){
        
        #function to pick up the most frequent label
        zuida<-0
        value.max<-c()
        geshu<-length(attr(table(x), "names"))
	
        for (b in 1:geshu){
            value<-attr(table(x), "names")[b]
        
            if (table(x)[[value]] > zuida) {
                zuida<-table(x)[[value]]
		value.max<-value
                } 
            }
	return (value.max)
        }

    return (apply(output, 1, major))

    }





combiVote<-function(data, known.index=1:nrow(data), 
                    k=ncol(data), N=50, npairs=100, 
                    epsilon=0.001, ...){
        
        #known.index is the index for instances with labels
	#check input
	if (is.matrix(data)) {
		data<-as.data.frame(data)
		#cat ("Ok!\n")
	} else  if (is.data.frame(data)){
		#cat ("good!\n")
	} else {
		stop ("Input must be either matrix or dataframe!")
		}	

	#discretize the data, train and test data together!
        pair<-c()
        triple<-c()
        disc.fea<-c()

        for (j in 1:N){
	    er<-sample((1:ncol(data))[-k], 2)
            pair<-rbind(pair, er)
            disc.fea<-cbind(disc.fea,Disc.kmeans(as.matrix(data[,er]), k=9, times=100))

            san<-sample((1:ncol(data))[-k], 3)
            triple<-rbind(triple, san)
	    disc.fea<-cbind(disc.fea, Disc.kmeans(as.matrix(data[,san]), k=25, times=100))
            }


        #partition data into train and test
        #check this!
	labels<-names(table(data[known.index,k]))
	train1_index<-c()
	for (i in 1:length(labels)){
		num<-nrow(data[data[,k]==labels[i],])
		index<-which(data[,k]==labels[i])
		train1_index<-c(train1_index,
			sample(index, round(num/2, 0)))	
		
		}


	###################check above!##################
	#boostrap to get a better estimation of conditional probs
        train<-sample(train1_index, 10*length(train1_index), replace=T)
       
        #make multiple classifier based on feature pairs.
        pred<-c()
        for (q in 1:npairs){
        
        #100
        two.col<-sample((1:ncol(disc.fea)), 2) 
        
        classify<-function (f1.val, f2.val){
            biao<-table(disc.fea[train,two.col], data[train,k])

            #focus on binary classfication first!
            labels<-names(table(data[known.index,k]))

            probs<-c() 
            max.p<-c()
	    good.i<-c()
            for (i in 1:length(labels)){
	        #disc.fea[m,two.col[1]], disc.fea[m, two.col[2]],
                somevalue<-biao[f1.val, f2.val, labels[i]]/sum(biao[disc.fea[,,labels[i]]])
                probs<-c(probs, somevalue)
	        if (max.p < somevalue) {
                    max.p<-somevalue
                    good.i<-i
                    }   
                }
        
            return (sample(labels, 1, prob=probs))
            }

           pred<-cbind(pred, apply(disc.fea[, two.col], 1, classify))
	   }
     return(pred)
    }




###################################################################
###Use lasso regression to boost multiple classfiers
###Aug 12, 2014
###Aug 17, 2014
##################################################################
# require (glmnet)
pred<-c("N","N","N","P")
pred<-rbind(pred, rep("P", 4))
pred<-rbind(pred, c("N","N","P","P"))

know.index<-1:3

real.labels<-c("N","P","P")

#######################
###Aug11, 2014, Live For yourself first!
###Use lasso regression for boost
#######################

report<-function (pred, known.index, real.labels, 
                 output=c("class", "prob")){
    #if match.arg(output)=="class"    
    library(glmnet)
        
    #assigh dummy variable to pred
    
    pred.mat<-matrix(0, nrow(pred), ncol(pred))
    
    pred.mat[pred=="posi"] <- 1
    pred.mat[pred=="nega"] <- -1 
   
    cv.glmmod <- cv.glmnet(x=pred.mat[known.index,], y=as.factor(real.labels),alpha=1, family="binomial")
    
    best_lambda <- cv.glmmod$lambda.min
    
    glmmod<-glmnet(x=pred.mat[known.index,], y=as.factor(real.labels),alpha=1,family='binomial')
    
    corrected.pred<-predict(glmmod, s=best_lambda, newx=pred[-known.index,])
    
    if (output=="class"){
        return (corrected.pred)
        }
 
    else if (output == "prob") {
        }

    }



#http://stats.stackexchange.com/questions/72251/an-example-lasso-regression-using-glmnet-for-binary-outcome
#http://blog.revolutionanalytics.com/2013/05/hastie-glmnet.html
#http://www-stat.wharton.upenn.edu/~stine/mich/lasso.R
