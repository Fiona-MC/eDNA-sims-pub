# from sim output, I want to be able to run a logistic model
# then I want to be able to send it through the rest of the pipeline for INLA output to analyze whether it "got things right"

library(stats)

thisDir <- "/space/s1/fiona_callahan/multiSim5/"
runs <- 1:1

for (run in runs) {
    # load sitetab
    sim_sitetab_sampled <- read.csv(file = paste0(thisDir, "randomRun", run, "/sim_sitetab_sampled.csv"), header = TRUE)
    
    # TODO this isn't done -- gotta store these and decide whether they're "right"
    # labID is location but it should not be treated as numeric in this setup -- maybe just include Lat+Long
    model_sp1 <- glm(Sp1 ~ 1+Cov1+Cov2+Cov3+Age+Lat+Long+Sp2+Sp3, family=binomial(link='logit'), data=sim_sitetab_sampled)
    summary(model_sp1)
    model_sp2 <- glm(Sp2 ~ 1+Cov1+Cov2+Cov3+Age+Lat+Long+Sp1+Sp3, family=binomial(link='logit'), data=sim_sitetab_sampled)
    summary(model_sp2)
    model_sp3 <- glm(Sp3 ~ 1+Cov1+Cov2+Cov3+Age+Lat+Long+Sp1+Sp2, family=binomial(link='logit'), data=sim_sitetab_sampled)
    summary(model_sp3)
}

