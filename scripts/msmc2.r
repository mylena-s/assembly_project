setwd("/home/mylequis/Documents/miscelaneos/psmc/deepvariant")

pdf("LG11.all_effective_population_size_new.pdf")
mu <- 3.5e-9
gen <- 3
Dat<-read.table("LG11.out.msmc2.final.txt", header=TRUE)
plot(Dat$left_time_boundary/mu*gen, (1/Dat$lambda)/(2*mu), log="x",ylim=c(0,300000),
     type="n", xlab="Years ago", ylab="effective population size")
lines(Dat$left_time_boundary/mu*gen, (1/Dat$lambda)/(2*mu), type="s", col="red")
legend("topright",legend=c("Ctuca"), col=c("red", lty=c(1,1)))
dev.off()