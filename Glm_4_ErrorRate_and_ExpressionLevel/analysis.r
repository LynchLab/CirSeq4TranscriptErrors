files<-list("FPKM_errors_Athaliana.txt","FPKM_errors_Atumefaciens.txt","FPKM_errors_Bsubtilis.txt","FPKM_errors_Ccrescentus.txt","FPKM_errors_Celegans.txt","FPKM_errors_Creinhardtii.txt","FPKM_errors_Dmelanogaster.txt","FPKM_errors_Dradiodurans.txt","FPKM_errors_Ecoli.txt","FPKM_errors_Hsapiens.txt","FPKM_errors_Hvolcanii.txt","FPKM_errors_Kradiotolerans.txt","FPKM_errors_Mflorum.txt","FPKM_errors_Mmusculus.txt","FPKM_errors_Msmegmatis.txt","FPKM_errors_Pcaudatum.txt","FPKM_errors_Pfluorescens.txt","FPKM_errors_Ptetraurelia.txt","FPKM_errors_Saureus.txt","FPKM_errors_Scerevisiae.txt","FPKM_errors_Senterica.txt")


sink("glm_out.txt", append = FALSE)

for (x in files) {
        data<-read.table(x,header=TRUE)
        data$CovlogFPKM<-data$Cov*log(data$FPKM)
        print(x)
        starting_vals = c(10^-7)
        model_noslope_glm<-glm(Num_errors ~ Cov +0,start = starting_vals, family = poisson(link = identity), data)
        starting_vals = c(10^-6, 10^-7)
        model_glm<-glm(Num_errors ~ Cov + CovlogFPKM +0,start = starting_vals, family = poisson(link = identity), data)
        print(summary(model_glm))
        print(anova(model_noslope_glm, model_glm, test = "Chisq"))
}

sink()
