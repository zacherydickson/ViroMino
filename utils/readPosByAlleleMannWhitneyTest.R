#!/usr/bin/env Rscript

maxPrecision = 7;
maxPhred=255;


lowerHalf <- function(x){
    ifelse(x > 0.5, 1 - x, x)
}

calcMannWhitneyP <- function(x,idx,jdx){
    suppressWarnings(wilcox.test(x[[idx]],x[[jdx]])$p.value)
}

phredScale <- function(p){
    (-10 * log10(p)) |> floor()
}

###MAIN

args <- commandArgs(trailingOnly = T)

if(length(args) < 1){
    stop("Requires path to file with readPositions")
}

file <- args[1]

df <- suppressWarnings(read.csv(file,sep="\t"))

posList <- split(df$Positions,df$AlleleID) |>
    lapply(strsplit,",") |>
    lapply("[[",1) |>
    lapply(as.numeric) |> 
    lapply(round,maxPrecision) |>
    lapply(lowerHalf)

if(length(posList) < 2 | is.null(posList$Ref)){
    cat(maxPhred)
} else {
    results <- sapply(1:(length(posList)-1),calcMannWhitneyP,x=posList,jdx=length(posList)) |> phredScale()
    paste(results,collapse=",") |> cat()
}


