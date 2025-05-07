#!/usr/bin/env Rscript

maxPrecision = 7;
maxPhred=255;
MU=0.25

lowerHalf <- function(x){
    ifelse(x > 0.5, 1 - x, x)
}

calcMannWhitneyP <- function(x,idx,jdx){
    suppressWarnings(wilcox.test(x[[idx]],x[[jdx]])$p.value)
}

calcWilcoxSignedRankSum <- function(x,mu=MU){
    if(!length(x)){
        return(fromPhred(maxPhred))
    }
    wilcox.test(x,alternative="less",mu=mu)$p.value |> suppressWarnings()
}

phredScale <- function(p){
    (-10 * log10(p)) |> floor()
}

fromPhred <- function(ph){
    10^(ph/-10)
}

###MAIN

args <- commandArgs(trailingOnly = T)

if(length(args) < 1){
    stop("Requires path to file with readPositions")
}

file <- args[1]

df <- suppressWarnings(read.csv(file,sep="\t"))

#We use the lower half because reads maybe reverse oriented in their mapping,
#therefore we can't differentiate between the termini, instead we look at the distance from a terminus

posList <- split(df$Positions,df$AlleleID) |>
    lapply(strsplit,",") |>
    lapply("[[",1) |>
    lapply(as.numeric) |> 
    lapply(round,maxPrecision) |>
    lapply(lowerHalf)

#The Ref allele might have no reads supporting it if there are multiple alts, 
#It might also have no reads if it a mapping artifact
#
if(is.null(posList$Ref)){
    posList$Ref <- 0
}

results <- sapply(posList,calcWilcoxSignedRankSum) |> phredScale()
#Reorder the results so the reference comes first
results <- results[c(length(results),1:(length(results)-1))]
#paste things together
paste(results,collapse=",") |> cat()

#There is an assumption here that (with the possible exception of the reference) all alleles have at least one
# supporting read therefor the output string will be in the order ref,alt1,...


