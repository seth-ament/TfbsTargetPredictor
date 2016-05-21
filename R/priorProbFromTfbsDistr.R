# fit a model to predict whether a gene is regulated by a TF, given the distribution of TF binding sites

options(stringsAsFactors=F)
require( caret )

shRNA = read.delim("Cusanovich2014_TF_siRNA.txt")

DEp = grep( "DE_p" , colnames(shRNA) )
DEp = shRNA[,DEp]
rownames(DEp) = shRNA$Symbol
colnames(DEp) = gsub("\\_(.*)" , "" , colnames(DEp) )
kdtfs = colnames(DEp)

binding = grep("binding",colnames(shRNA))
binding = shRNA[ , binding ]
bind.tfs = gsub("\\_(.*)" , "" , colnames(binding) )
colnames(binding) = bind.tfs

logp = -log10(DEp)
tfs = intersect( colnames(logp) , colnames(binding) )

de.int = as.matrix( logp[ , tfs ] )
bind.int = as.matrix( binding[ , tfs ] )

a = as.vector( de.int )
b = as.vector( bind.int )
cor.test( a , b )

t = table( a > 3 , b >= 1 )


shRNA.genes = rownames(logp)

# modeled TFs
probdistr.files = Sys.glob("TfbsProbDistr/*")
modeled.tfs = gsub("TfbsProbDistr/","",probdistr.files)
modeled.tfs = gsub(".RData","",modeled.tfs)

tflist = intersect( kdtfs , modeled.tfs )

  tf = tflist[1]
  load( paste("TfbsProbDistr/",tf,".RData" , sep="" ) )
  tfbs = outp[ duplicated(outp$gene_name) == F , -c(1:3) ]
  de = data.frame( rownames(logp) , logp[,tf] )
  colnames(de)[2] = "logp"
  merged = merge( tfbs , de ,by.x = 1 , by.y = 1 )
  probdistr = cbind( tf , merged )
  colnames(probdistr)[1] = "tf"
for( tf in tflist[-1] ) {
  load( paste("TfbsProbDistr/",tf,".RData" , sep="" ) )
  tfbs = outp[ duplicated(outp$gene_name) == F , -c(1:3) ]

  de = data.frame( rownames(logp) , logp[,tf] )
  colnames(de)[2] = "logp"
  merged = merge( tfbs , de ,by.x = 1 , by.y = 1 )
  probdistr.tf = cbind( tf , merged )
  colnames(probdistr.tf)[1] = "tf"
  probdistr = rbind( probdistr , probdistr.tf )
}

x = rowSums( probdistr[,5:8] )
y = probdistr$logp > 3

tfs = intersect( colnames(de.int) , unique(probdistr$tf ))
tf = "E2F1"
a = rowSums( probdistr[ , 4:9 ] )
a = probdistr$gene_name[ a >= 2 & probdistr$tf == tf ]
b = shRNA$Symbol[ bind.int[,tf] >= 1 ]
g = unique( probdistr$gene_name )
t =table( g %in% a  , g %in% b ) 
fisher.test( t )
t

ta


fisher.test( table( g %in% a  , g %in% b ) )
 

t = table( rowSums(probdistr[,5:8]) > 1 , probdistr$logp > 3 ) 
fisher.test( t )

t.test( x ~ y )

merged = merge( probdistr , logp , by.x = 4 , by.y = 0 )

genes0 = merged[,1]

nas = is.na( merged[,tf] )

x = as.matrix( merged[ nas == F ,  6:11 ])
y = merged[ nas == F , tf ]

for( i in 1:6 ) {
  a = x[ y > 3 , i ] > 1
  b = x[ y < 1 , i ] > 1
  table( a , b )
  test = t.test( x[ y > 3 , i ] , x[ y < 1 , i ] ) 
  cat( test$p.value , "\n" )
}

x.tot = rowSums( x[2:4] )
t.test( x.tot ~ factor(y>3) )


ends = y > 3 | y < 1


inTrain = createDataPartition( y , p = 0.75 , list = F )
x.train = x[ inTrain , ]
y.train = y[ inTrain ]
x.test = x[ -inTrain , ]
y.test = y[ -inTrain ]

cor.test( y x[,3] )





trControl = trainControl( method = "cv" )
grid = expand.grid( n.trees = c(50,100,250,500) ,
		interaction.depth = c(1,3,5) ,
		shrinkage = 0.1 ,
		n.minobsinnode = 10 )
fit = train( x = x.train , y = factor( y.train > 2 ) ,
	trControl = trControl , metric = "Kappa" ,
	tuneGrid = grid , method = "gbm" )
pred.train = predict( fit , newdata = x.train , n.trees = 500 , type = "prob" )
pred = predict( fit , newdata = x.test , n.trees = 500 , type = "prob" )

require( glmnet )
fit2 = cv.glmnet( x = x.train , y = factor( y.train > 3 ) , 
	family = "binomial" , alpha = 0 )
pred = predict( fit2 , newx = x.train , s = fit2$lambda.1se )[,1]


table( y.train > 2 , pred.train[,2] > 0.5 )










