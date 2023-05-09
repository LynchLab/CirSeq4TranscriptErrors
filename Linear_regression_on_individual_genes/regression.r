files_individuals<-list("FPKM_ErrorRates_Athaliana.txt","FPKM_ErrorRates_Atumefaciens.txt","FPKM_ErrorRates_Bsubtilis.txt","FPKM_ErrorRates_Ccrescentus.txt","FPKM_ErrorRates_Celegans.txt","FPKM_ErrorRates_Creinhardtii.txt","FPKM_ErrorRates_Dmelanogaster.txt","FPKM_ErrorRates_Dradiodurans.txt","FPKM_ErrorRates_Ecoli.txt","FPKM_ErrorRates_Hsapiens.txt","FPKM_ErrorRates_Hvolcanii.txt","FPKM_ErrorRates_Kradiotolerans.txt","FPKM_ErrorRates_Mflorum.txt","FPKM_ErrorRates_Mmusculus.txt","FPKM_ErrorRates_Msmegmatis.txt","FPKM_ErrorRates_Pcaudatum.txt","FPKM_ErrorRates_Pfluorescens.txt","FPKM_ErrorRates_Ptetraurelia.txt","FPKM_ErrorRates_Saureus.txt","FPKM_ErrorRates_Scerevisiae.txt","FPKM_ErrorRates_Senterica.txt")


sink("regression_out.txt", append = FALSE)


for (x in files_individuals) {
	data<-read.table(x,header=TRUE)
	model<-lm(formula = ErrorRate ~ log2(FPKM),data)
	print(x)
	print(summary(model))
}

sink()
