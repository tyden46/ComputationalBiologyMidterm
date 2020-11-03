setwd("C:/Users/tyson/OneDrive/Desktop/Fall 2020 Classes/Computational Biology Midterm/Problem 2")
library(genbankr)
library(Biostrings)
library(ggplot2)
id=GBAccession("NC_045512.2") #Accession number for SARS-CoV-2 Ref
gb=readGenBank(id)
#Read the genome from the genbank file
fullRefGenome=toString(gb@sequence$`Severe acute respiratory syndrome coronavirus2`)

id=GBAccession("MW192772.1") #Accession number for SARS-CoV-2 Ref
gb=readGenBank(id)

#Read the genome from the genbank file
fullAusGenome=toString(gb@sequence$`Severe acute respiratory syndrome coronavirus2`)

#Playing with Mismatch Scores
mismatchScores=-1*(0.5 * 1.5^(0:20))
alignmentScores=c()
for(x in mismatchScores){
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = x, baseOnly = TRUE)
  alignmentScore=pairwiseAlignment(fullAusGenome, fullRefGenome, substitutionMatrix=mat, scoreOnly=TRUE)
  alignmentScores=append(alignmentScores,alignmentScore)
}

mytable=(as.data.frame(cbind(log10(-mismatchScores),alignmentScores)))
colnames(mytable)=c("logOfTheMismatchScore")
png(filename="MismatchScore.png",width=500,height=500)
ggplot(mytable, aes(x=logOfTheMismatchScore, y=alignmentScores)) +
  geom_point() +
  ggtitle(paste("As the magnitude of the mismatch", "score increases the",
                 "Alignment score decreases", "until it hits a minimum",
                 sep='\n'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
#As the magnitude of the mismatch score increases the 
#Alignment score decreases until it hits a minimum

#Playing with Gap Scores
gapScores=0.5 * 1.5^(0:20)
alignmentScores=c()
for(x in gapScores){
  alignmentScore=pairwiseAlignment(fullAusGenome, fullRefGenome, gapExtension=x, scoreOnly=TRUE)
  alignmentScores=append(alignmentScores,alignmentScore)
}

mytable=(as.data.frame(cbind(log10(gapScores),alignmentScores)))
colnames(mytable)=c("logOfTheGapScore","Alignment Score")

png(filename="GapScore.png",width=500,height=500)
ggplot(mytable, aes(x=logOfTheGapScore, y=`Alignment Score`)) +
  geom_point()+
  ggtitle(paste("As the magnitude of the gap", "score increases the",
                "Alignment score decreases", "and does not hit a minimum",
                sep='\n'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()