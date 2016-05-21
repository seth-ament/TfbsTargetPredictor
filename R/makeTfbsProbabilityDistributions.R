# sum the tfbs probabilities around each TSS

require( TReNA )
require( GenomicRanges )

# load the 
tfsWithModels = Sys.glob("/proj/price1/sament/lymphoblast_trn/FittedTFBSProbabilities/*")
tfsWithModels = gsub(".RData","",tfsWithFootprints)
tfsWithModels = gsub("/proj/price1/sament/lymphoblast_trn/FittedTFBSProbabilities/","",tfsWithFootprints)

tfsWithFootprints = Sys.glob("/proj/price1/sament/lymphoblast_trn/AllFeatureMatrix/*")
tfsWithFootprints = gsub("/proj/price1/sament/lymphoblast_trn/AllFeatureMatrix/","",tfsWithFootprints)
tfsWithFootprints = gsub(".RData","",tfsWithFootprints)

tfsWithChIP = Sys.glob("/proj/price1/sament/lymphoblast_trn/known_tfbs/hg38/*")

tfsWithChIP = gsub("/proj/price1/sament/lymphoblast_trn/known_tfbs/hg38/","",tfsWithChIP)

tfsWithChIP = gsub("\\_(.*)" , "" , tfsWithChIP )


tflist = intersect( tfsWithChIP , intersect( tfsWithModels , tfsWithFootprints ))

# get transcript models and TSSs

# open database conn
genome.db.uri <- "postgres://whovian/hg38"
project.db.uri <-  "postgres://whovian/lymphoblast"
fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
# get transcripts from gtf
query <-
paste( "select gene_name, transcript_name, chr, start, endpos, strand from gtf where" ,
        "gene_biotype='protein_coding' and moleculetype='transcript'" , sep=" " )

transcripts = dbGetQuery( fp@genome.db , query )
# close db conn
closeDatabaseConnections(fp)

# function to get each transcript's TSS
get_tss <-
function( t ) {
      chrom <- transcripts$chr[t]
      start.orig <- transcripts$start[t]
      end.orig   <- transcripts$endpos[t]
      strand     <- transcripts$strand[t]

      if(strand == "-"){ # reverse (minus) strand.  TSS is at "end" position
         tss <- end.orig
         }
      else{ #  forward (plus) strand.  TSS is at "start" position
        tss <- start.orig
        }
     return( tss )
}
# apply get_tss to all transcripts
tss = sapply( 1:nrow(transcripts) , get_tss )
# assemble a bed file for the TSSs
tss.loc = unique( data.frame(
        chr = transcripts$chr ,
        start = tss , end = tss ,
        gene_name = transcripts$gene_name ,
	transcript_name = transcripts$transcript_name ) )
tss.loc = tss.loc[ duplicated( tss.loc[,1:3] ) == F , ]
# GRanges obj
tss.gr = makeGRangesFromDataFrame( tss.loc )

require(foreach)
require(doParallel)
registerDoParallel( cores = 8 )

# calcualte the precision and recall for each TF with ChIP-seq data

prec.recall =
foreach( tf=tflist ) %dopar% {

cat("working on",tf,"\n")

prob_file = paste("FittedTFBSProbabilities/",tf,".RData",sep="")
load( prob_file )
probs = bed[,-c(1:7)]

chip = Sys.glob( paste(
  "/proj/price1/sament/lymphoblast_trn/known_tfbs/hg38/",tf,"*",sep=""))
chip = read.table( chip )
colnames(chip)[1:3] = c("chr","start","end")
chip = makeGRangesFromDataFrame(chip)

o3 = countOverlaps( tss.gr , chip , maxgap = 10000 )
predictor.perc = seq(0.3,0.9,0.1) #0.4,0.5,0.6,0.7,0.8,0.9)
prob.thresh = seq(0,0.9,0.1)
m = length(predictor.perc)
n = length(prob.thresh)
precision.10kb = matrix( NA , m ,n)
recall.10kb = matrix( NA , m ,n)
precision.tfbs = matrix( NA , m ,n)
recall.tfbs = matrix( NA , m ,n)
colnames(precision.10kb) = colnames(recall.10kb) = 
  colnames(precision.tfbs) = colnames(recall.tfbs) = paste("Pr>",prob.thresh,sep="")
rownames(precision.10kb) = rownames(recall.10kb) = 
  rownames(precision.tfbs) = rownames(recall.tfbs) = paste("k>",predictor.perc,sep="")
for( i in 1:m ) {
 cat(i,":")
 prob = apply( probs , 1 , quantile , probs = predictor.perc[i] )
 for( j in 1:n ) {
  cat(j)
  x = which( prob > prob.thresh[j] )
  if( length(x) == 0 ) next
  fp = makeGRangesFromDataFrame( bed[ x , ] )
  o1 = countOverlaps( chip , fp , maxgap = 100 )
  o2 = countOverlaps( fp , chip , maxgap = 100 )
  precision.tfbs[i,j] = length(which(o2>0)) / length(o2)
  recall.tfbs[i,j] = length(which( o1>0)) / length(o1)
  o4 = countOverlaps( tss.gr , fp , maxgap = 10000 )
  if( all( o4 == 0 ) ) next
  t = table( o3 > 0 , o4 > 0 )
  recall.10kb[i,j] = t[2,2] / sum(t[2,])
  precision.10kb[i,j] = t[2,2] / sum(t[,2])
 }
 cat("\n")
}

return( list( 	precision.tfbs = precision.tfbs ,
		recall.tfbs = recall.tfbs ,
		precision.10kb = precision.10kb ,
		recall.10kb = recall.10kb ) )
}

names( prec.recall ) = tflist
save( prec.recall , file="precision_and_recall_of_tfbs_predictions_for_57_tfs.RData" )

# find an optimal balance between precision and recall
# across the 57 TFs
# formula: minimize the sum of the dropoff in precision and the dropoff in recall
# use a modifier weight (called 'ratio') to weight precision vs. recall
# ratio * ( max(precision) / precision ) + ( max(recall) / recall )

ratio = 1.5
pr = 
foreach( tf = tflist ) %dopar% {
  p.tfbs = prec.recall[[tf]]$precision.tfbs
  r.tfbs = prec.recall[[tf]]$recall.tfbs
  pr.tfbs = (ratio*(max(p.tfbs,na.rm=T)/p.tfbs)) + (max(r.tfbs,na.rm=T)/r.tfbs)
  p.10kb = prec.recall[[tf]]$precision.10kb
  r.10kb = prec.recall[[tf]]$recall.10kb
  pr.10kb = (1.5*(max(p.10kb,na.rm=T)-p.10kb)) + (max(r.10kb,na.rm=T)-r.10kb)
  return( list( pr.tfbs = pr.tfbs ,
		pr.10kb = pr.10kb ) )
}
names(pr) = tflist

# rank the values in the grid based on the combined 'pr' scores
rank.pr.tfbs =
foreach( tf = tflist ) %dopar% {
  order(pr[[tf]]$pr.tfbs)
}
rank.pr.tfbs = as.data.frame( rank.pr.tfbs )
colnames(rank.pr.tfbs) = tflist
# identify the optimal values based on the grid ranks for all 57 TFs
matrix( apply( rank.pr.tfbs , 1 , median ) , 7 , 10 )

rank.pr.10kb =
foreach( tf = tflist ) %dopar% {
  order(pr[[tf]]$pr.10kb)
}
rank.pr.10kb = as.data.frame( rank.pr.10kb )
colnames(rank.pr.10kb) = tflist
matrix( apply( rank.pr.10kb , 1 , median ) , 7 , 10 )


# best values for predicting each TFBS:
# use the 3rd highest probability estimate across the 'ensemle' of models
# consider as TFBSs those that have a probability > 0.3

# best values for counting the number of TFBSs in the 10kb window around each TSS:
# use the 2nd highest probability estimate across the ensemble of models
# consider as TFBSs those that have a probability > 0.2







# code for calculating the distribution of TFBS probabilities
# around each TSS. The method used in this version did not 
# perform well.

left = c(-10^(5:3),0,10^(3:4))
right = c(-10^(4:3),0,10^(3:5))

tfbs.prob.distribution <- 
lapply( 1:length(left) , function(r) {
  # cat("working on" , r , "\n" )
  regions.left = tss.loc$start + left[r]
  regions.left[ regions.left < 0 ] = 0
  regions.right = tss.loc$start + right[r]
  regions.right[ regions.right < 0 ] = 1
  chrom = tss.loc$chr
  df = data.frame( chrom , start = regions.left , end = regions.right )
  gr = makeGRangesFromDataFrame( df )
  hits = as.matrix( findOverlaps( gr , footprints ))
  tx.ind = unique( hits[,1] )
  prob.sum = rep( 0 , nrow(tss.loc) )
  prob.sum[ tx.ind ] <-
  sapply( tx.ind , function(i) {
    fp.ind = hits[ hits[,1] == i , 2 ]
    probs = prob.tfbs[ fp.ind ]
    return( sum(probs) )
  }
  )
  return( prob.sum )
} 
) # end internal loop (loop over regions)

tfbs.prob.distribution = as.data.frame( tfbs.prob.distribution )
colnames(tfbs.prob.distribution) = c( 
	"minus1e5_1e4" ,
	"minus1e4_1e3" ,
	"minus1e3_0" ,
	"plus0_1e3" ,
	"plus1e3_1e4" ,
	"plus1e4_1e5" )


outp = cbind( tss.loc , tfbs.prob.distribution )
save( outp , file=paste("TfbsProbDistr/",tf,".RData",sep="") )
rm( prob_file )
} # end main loop (looping over tfs)



