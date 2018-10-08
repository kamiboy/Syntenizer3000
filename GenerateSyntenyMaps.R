if(!require(MASS)){
  install.packages("MASS")
  library(MASS)
}
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1)
{
  stop("When running this script please provide path of folder containing map data (the output folder of Syntenizer3000) as an argument.", call.=FALSE)
} else
{
  setwd(args[1])
  if (file.exists("gene_groups.csv"))
  {
    groups<-read.csv(file="gene_groups.csv", sep=":", header = F)[,1]
    print("Rendering synteny maps, please wait...")
    progress <- -1
    for (group in 1:length(groups)){
      map = read.csv(file = paste(groups[group], "_synteny_chart.csv", sep=""), dec=",", sep=";", header=T, check.names = F)
      strains = colnames(map)
      
      pdf(paste(groups[group], "_synteny_chart.pdf", sep=""), width=15, height=15)
      par(las=2, cex=0.5, lwd=0.25)
      parcoord(map, lty = 1, main=groups[group], col= c("#00000000", "#00000000", "black", sample(rainbow(nrow(map)-3))))
      dev.off()
      
      if (progress != as.integer( (group*10)/length(groups)))
      {
        progress = as.integer((group*10)/length(groups))
        print(paste(progress*10, "%", sep="") )
      }
    }
    print("Rendering complete.")
  }
  else
  {
    stop("Could not locate gene_groups.csv", call.=FALSE)
  }
}

