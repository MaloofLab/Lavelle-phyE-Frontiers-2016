## reducing data to what is really used in the paper

### Figure 2 data

physunshade <- read.csv("physunshade.csv")
summary(physunshade)
levels(physunshade$genotype)

physunshade.small <- droplevels(physunshade[grepl("(Ler)|(phyBD$)|(moneymaker)|(phyB1B2)",physunshade$genotype),])

physunshade.small

physunshade.small$species <- factor(physunshade.small$species,labels=c("Arabidopsis","S.lyc"))

summary(physunshade.small)

write.csv(physunshade.small,"Fig2AB_hyps.csv")

phymutantchamber <- read.csv("phymutantchamber.csv")

summary(phymutantchamber)

levels(phymutantchamber$genotype)

phymutantchamber.small <- droplevels(phymutantchamber[
  phymutantchamber$week==5 & 
    grepl("(Moneymaker)|(phyB2/phyB1)",phymutantchamber$genotype),
                                           c(1:8,14)])

summary(phymutantchamber.small)

write.csv(phymutantchamber.small,"Fig2C_tomato_int.csv")

phyintpetbothreps <- read.csv("phyintpetbothreps.csv")

phyintpetbothreps <- phyintpetbothreps[!grepl("3",phyintpetbothreps$genotype),]
phyintpetbothreps <- droplevels(phyintpetbothreps)