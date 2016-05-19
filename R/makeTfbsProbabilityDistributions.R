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

tflist = intersect( tfsWithModels , tfsWithFootprints )

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

left = c(-10^(5:3),0,10^(3:4))
right = c(-10^(4:3),0,10^(3:5))

require(foreach)
require(doParallel)
registerDoParallel( cores = 12 )

# get the footprints and their probabilities for one TF

foreach( tf=tflist ) %dopar% {
cat("working on",tf,"\n")

prob_file = paste("FittedTFBSProbabilities/",tf,".RData",sep="")
load( prob_file )
prob.tfbs = bed$median.probs
footprints = makeGRangesFromDataFrame( bed ) # , keep.extra.columns = T )

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



