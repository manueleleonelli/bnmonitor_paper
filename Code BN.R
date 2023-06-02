library(bnlearn)
library(qgraph)
library(bnmonitor)
library(dplyr)

data <- readRDS(file = "data.rds")
setwd("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis")
diabetes <- read.csv("diabetes.csv")
colnames(diabetes) <- c("PREG","GLUC","PRES","TRIC","INS","MASS","PED","AGE","DIAB")
diabetes <- data.frame(apply(diabetes, 2, as.numeric))
diabetes[,-9] <- bnlearn::discretize(diabetes[,-9], "quantile", 2)
diabetes$PREG <- factor(diabetes$PREG, labels = c("low","high"))
diabetes$GLUC <- factor(diabetes$GLUC,labels = c("low","high"))
diabetes$PRES <- factor(diabetes$PRES, labels = c("low","high"))
diabetes$TRIC <- factor(diabetes$TRIC,labels = c("low","high"))
diabetes$INS <- factor(diabetes$INS, labels = c("low","high"))
diabetes$MASS <- factor(diabetes$MASS,labels = c("low","high"))
diabetes$PED <- factor(diabetes$PED, labels = c("low","high"))
diabetes$AGE <- factor(diabetes$AGE,labels = c("low","high"))
diabetes$DIAB <- factor(diabetes$DIAB,levels =c(0,1),labels=c("neg","pos"))


data <- diabetes
data <- data[,c("AGE","PREG","MASS","GLUC","PRES","TRIC","INS","PED","DIAB")]
#diabetes <- read.csv("diabetes.csv")
#data<- read.delim("C:/Users/manuele.leonelli/Downloads/customer_survey.txt")
## Create train and test data
set.seed(2001)
indices <- sample(nrow(data), 450)
train <- data[indices,]
#test <- test[order(test$REV),]
test <- data[-indices,]


## BLACKLIST
#bl <- cbind(rep("GROWTH",14),colnames(data)[1:14])

bl <- rbind(expand.grid("DIAB",colnames(data)[-9]),expand.grid("GLUC",colnames(data)[1:3]),expand.grid("PRES",colnames(data)[1:3]),expand.grid("TRIC",colnames(data)[1:3]),expand.grid("INS",colnames(data)[1:3]),expand.grid("PED",colnames(data)[1:3]), expand.grid("MASS",colnames(data)[1:2]), expand.grid("PREG","AGE") )

layout <- matrix(c(2,0,0,0,1,1,1,2,0.5,0.5,0.5,1.5,0,2,1.5,1.5,2,2),byrow=T,nrow=9)

### LEARN DAG TABU-BIC

boot.tabu.bic <- boot.strength(train, R = 1000, m = nrow(train), algorithm = "tabu", algorithm.args = list(score = "bic", blacklist = bl), cpdag = F, debug = FALSE)
arc.set <- boot.tabu.bic [(boot.tabu.bic$strength > 0.7) & (boot.tabu.bic$direction >= 0.5), ]
dag.tabu.bic <- bnlearn::empty.graph(colnames(train))
bnlearn::arcs(dag.tabu.bic) <- arc.set[,1:2]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/bic.png", units="cm", width=12, height=6, res=300)
qgraph(dag.tabu.bic,layout=layout, edge.width= 4, vsize = 10)
dev.off()

### LEARN DAG TABU-AIC
boot.tabu.aic <- boot.strength(train, R = 1000, m = nrow(train), algorithm = "tabu", algorithm.args = list(score = "aic", blacklist = bl), cpdag = F, debug = FALSE)
arc.set <- boot.tabu.aic [(boot.tabu.aic$strength > 0.7) & (boot.tabu.aic$direction >= 0.5), ]
dag.tabu.aic <- bnlearn::empty.graph(colnames(train))
bnlearn::arcs(dag.tabu.aic) <- arc.set[,1:2]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/aic.png", units="cm", width=8, height=5, res=300)
qgraph(dag.tabu.aic,layout=layout,edge.width= 4, vsize = 10)
dev.off()

### LEARN DAG TABU-FNML
boot.tabu.fnml <- boot.strength(train, R = 1000, m = nrow(train), algorithm = "tabu", algorithm.args = list(score = "fnml", blacklist = bl), cpdag = F, debug = FALSE)
arc.set <- boot.tabu.fnml [(boot.tabu.fnml$strength > 0.7) & (boot.tabu.fnml$direction >= 0.5), ]
dag.tabu.fnml <- bnlearn::empty.graph(colnames(train))
bnlearn::arcs(dag.tabu.fnml) <- arc.set[,1:2]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/fnml.png", units="cm", width=8, height=5, res=300)
qgraph(dag.tabu.fnml,layout=layout,edge.width= 4, vsize = 10)
dev.off()

### LEARN DAG TABU-K2
boot.tabu.k2 <- boot.strength(train, R = 1000, m = nrow(train), algorithm = "tabu", algorithm.args = list(score = "k2", blacklist = bl), cpdag = F, debug = FALSE)
arc.set <- boot.tabu.k2 [(boot.tabu.k2$strength > 0.7) & (boot.tabu.k2$direction >= 0.5), ]
dag.tabu.k2 <- bnlearn::empty.graph(colnames(train))
bnlearn::arcs(dag.tabu.k2) <- arc.set[,1:2]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/k2.png", units="cm", width=8, height=5, res=300)
qgraph(dag.tabu.k2,layout=layout,edge.width= 4, vsize = 10)
dev.off()

### LEARN DAG TABU-BDLA
boot.tabu.bdla <- boot.strength(train, R = 1000, m = nrow(train), algorithm = "tabu", algorithm.args = list(score = "bdla", blacklist = bl), cpdag = F, debug = FALSE)
arc.set <- boot.tabu.bdla [(boot.tabu.bdla$strength > 0.7) & (boot.tabu.bdla$direction >= 0.5), ]
dag.tabu.bdla <- bnlearn::empty.graph(colnames(train))
bnlearn::arcs(dag.tabu.bdla) <- arc.set[,1:2]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/bdla.png", units="cm", width=8, height=5, res=300)
qgraph(dag.tabu.bdla,layout=layout,edge.width= 4, vsize = 10)
dev.off()

### LEARN DAG GS
boot.gs <- boot.strength(train, R = 1000, m = nrow(train), algorithm = "gs", algorithm.args = list(blacklist = bl), cpdag = F, debug = FALSE)
arc.set <- boot.gs [(boot.gs$strength > 0.7) & (boot.gs$direction >= 0.5), ]
dag.gs <- bnlearn::empty.graph(colnames(train))
bnlearn::arcs(dag.gs) <- arc.set[,1:2]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/gs.png", units="cm", width=8, height=5, res=300)
qgraph(dag.gs,layout=layout,edge.width= 4, vsize = 10)
dev.off()

### LEARN DAG PC
boot.iamb <- boot.strength(train, R = 1000, m = nrow(train), algorithm = "pc.stable", algorithm.args = list(blacklist = bl), cpdag = F, debug = FALSE)
arc.set <- boot.iamb [(boot.iamb$strength > 0.7) & (boot.iamb$direction >= 0.5), ]
dag.iamb <- bnlearn::empty.graph(colnames(train))
bnlearn::arcs(dag.iamb) <- arc.set[,1:2]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/pc.png", units="cm", width=8, height=5, res=300)
qgraph(dag.iamb,layout=layout,edge.width= 4, vsize = 10)
dev.off()

### GLOBAL AND NODE MONITOR
dags <- list(dag.tabu.bic,dag.tabu.aic,dag.tabu.fnml,dag.tabu.k2,dag.tabu.bdla,dag.gs,dag.iamb)
sapply(dags,global_monitor,df = test)
lapply(dags,node_monitor,df = test)


## SEQUENTIAL MONITORS
test1<- test[order(diabetes$Insulin),]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/cond_bic.png", units="cm", width=10, height=5, res=300)
plot(seq_cond_monitor(dag.tabu.bic,test1,"GLUC"))
dev.off()
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/cond_gs.png", units="cm", width=10, height=5, res=300)
plot(seq_cond_monitor(dag.gs,test1,"GLUC"))
dev.off()
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/marg_bic.png", units="cm", width=10, height=5, res=300)
plot(seq_marg_monitor(dag.tabu.bic,test1,"GLUC"))
dev.off()
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/marg_gs.png", units="cm", width=10, height=5, res=300)
plot(seq_marg_monitor(dag.gs,test1,"GLUC"))
dev.off()

test1<- test[order(diabetes$SkinThickness),]
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/cond_bic1.png", units="cm", width=10, height=5, res=300)
plot(seq_cond_monitor(dag.tabu.bic,test1,"DIAB"))
dev.off()
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/cond_gs1.png", units="cm", width=10, height=5, res=300)
plot(seq_cond_monitor(dag.gs,test1,"DIAB"))
dev.off()
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/marg_bic1.png", units="cm", width=10, height=5, res=300)
plot(seq_marg_monitor(dag.tabu.bic,test1,"DIAB"))
dev.off()
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/marg_gs1.png", units="cm", width=10, height=5, res=300)
plot(seq_marg_monitor(dag.gs,test1,"DIAB"))
dev.off()

png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/high.png", units="cm", width=10, height=6, res=300)
plot(seq_pa_ch_monitor(dag.gs,test1,"DIAB","GLUC","high"))
dev.off()
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/low.png", units="cm", width=10, height=6, res=300)
plot(seq_pa_ch_monitor(dag.gs,test1,"PRES","GLUC","low"))
dev.off()

## INFLUENCE
influence <- influential_obs(dag.tabu.bic,test)
png("C:/Users/manuele.leonelli/Dropbox/bnmonitor/Analysis/img/influence.png", units="cm", width=11, height=5, res=600)
plot(influence)
dev.off()


### SENSITIVITY
bn <- bnlearn::bn.fit(dag.tabu.bic, diabetes)
sens_g <- sensitivity (bn, interest_node = "DIAB" , interest_node_value = "pos" , node = "GLUC" ,
                          value_node = "high" , value_parents = "high" ,new_value = "all" )
plot(sens_g)

sens_m <- sensitivity( bn , interest_node = "PREG" ,
                          interest_node_value = "high" ,
                          evidence_nodes = "DIAB" , evidence_states = "pos" ,
                          node = "MASS" , value_node = "high" ,
                          value_parents = NULL , new_value = "all" )
plot(sens_m )

cd_g <- CD( bn , node = "GLUC" , value_node = "high" ,
               value_parents = "high" , new_value = "all" )
plot (cd_g)
cd_m <- CD( bn , node = "MASS" , value_node = "high" ,
               value_parents = NULL , new_value = "all" )
plot (cd_m)

bnlearn::cpquery( bn , event = ( DIAB == "pos" ) ,
                     evidence = ( PRES == "high" ))

sensquery( bn , interest_node = "DIAB" ,
            interest_node_value = "pos" , new_value = 0.4 ,
            evidence_nodes = "PRES" , evidence_states = "high" )
