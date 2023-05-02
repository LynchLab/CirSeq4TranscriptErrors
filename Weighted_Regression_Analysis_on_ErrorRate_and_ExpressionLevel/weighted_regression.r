files_bins<-list("ExressionLevel_ErrorRate_20bins_Athaliana.txt","ExressionLevel_ErrorRate_20bins_Atumefaciens.txt","ExressionLevel_ErrorRate_20bins_Bsubtilis.txt","ExressionLevel_ErrorRate_20bins_Ccrescentus.txt","ExressionLevel_ErrorRate_20bins_Celegans.txt","ExressionLevel_ErrorRate_20bins_Creinhardtii.txt","ExressionLevel_ErrorRate_20bins_Dmelanogaster.txt","ExressionLevel_ErrorRate_20bins_Dradiodurans.txt","ExressionLevel_ErrorRate_20bins_Ecoli.txt","ExressionLevel_ErrorRate_20bins_Hsapiens.txt","ExressionLevel_ErrorRate_20bins_Hvolcanii.txt","ExressionLevel_ErrorRate_20bins_Kradiotolerans.txt","ExressionLevel_ErrorRate_20bins_Mflorum.txt","ExressionLevel_ErrorRate_20bins_Mmusculus.txt","ExressionLevel_ErrorRate_20bins_Msmegmatis.txt","ExressionLevel_ErrorRate_20bins_Pcaudatum.txt","ExressionLevel_ErrorRate_20bins_Pfluorescens.txt","ExressionLevel_ErrorRate_20bins_Ptetraurelia.txt","ExressionLevel_ErrorRate_20bins_Saureus.txt","ExressionLevel_ErrorRate_20bins_Scerevisiae.txt","ExressionLevel_ErrorRate_20bins_Senterica.txt")


sink("regression_out.txt", append = FALSE)


for (x in files_bins) {
	data<-read.table(x,header=TRUE)
	model<-lm(formula = ErrorRate ~ log2(ExpressionLevel),weights=Weights,data)
	print(x)
	print(summary(model))
}

sink()
