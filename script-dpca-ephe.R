##### Enlever les objets existants ##### 
rm(list=ls())

##### Installer les libraries R nécessaires pour le script ##### 
library(adegenet)
library(vcfR)
library(ggplot2)

################ 1 - Données génomiques ###################

### Importer les données vcf dans l'environnement de travail
SEBASTE_vcf <- read.vcfR('100snps_861ind_sebastes.recode.vcf')

### Transformer les données object vcf en object de type genind
SEBASTE_genind <- vcfR2genind(SEBASTE_vcf)

### Vérifier le nombre d'échantillons
SEBASTE_genind

### Importer les informations sur chaque échantillon (localité, latitude, longitude)
pop <- read.table('population_map_sebastes.txt',header=TRUE)


################# 2 - Analyse de regroupement ###################

# Déterminer le nombre de groupes potentiels via le critère d'information bayésien (en anglais bayesian information criterion ; en abrégé BIC) 
# Le BIC est un critère d'information dérivé du critère d'information d'Akaike proposé par Gideon Schwarz en 1978.
# Le modèle qui sera sélectionné est celui qui minimise le critère BIC.
grp <- find.clusters(SEBASTE_genind, max.n.clust=10)

# Choisissez le nombre d'axes à garder à partir duquel la courbe suit un plateau. Cependant ne pas sélectionner un nombre d'axes > N/3
# Choisir le nombre optimal de groupes génétiques (la plus petite valeur de BIC) quand le BIC décline de manière très brutale

### Observer la variation de valeur d'AIC en fonction du nombre de groupes
grp$Kstat

### Observer le nombre d'individus par groupe génétique
grp$size

# Regarder quels individus appartiennent à quel groupe génétique et est ce que ça fait du sens?
individuals <- head(grp$grp, 861) 
individuals

# Sauvegarder l'information pour chaque individus
data_kmeans <- data.frame(grp$grp)
data_kmeans 

# Exporter les données dans un fichier 
write.table(individuals, "Individuals_clusters_BIC.txt", quote=F)

################# 3 - DACP sans a prior ###################

### Faire tourner l'analyse de DPCA sans a priori, en prenant le regroupement trouvé par l'analyse de BIC
dapc_noprior <-dapc(SEBASTE_genind, grp$grp)

### Visualiser rapidemment les résultats de la DAPC
scatter(dapc_noprior, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

### Regarder la contribution de chaque SNP à la discrimination des groupes génétiques trouvés
names(dapc_noprior)
head(dapc_noprior$var.contr)
load_dapc_noprior <- as.data.frame(dapc_noprior$var.contr)

### Analyser le pourcentage de variation génétique observée à chacun des axes
pdf("Percent_dapc_all_loci.pdf")
percent= dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", names.arg=round(percent,2))
dev.off()
percent

################## Visualiser les résultas de la DAPC avec ggplot #####################

### Write the results of the DPCA and modify these results on excel in order to get a file such as "PC_results_adegenet.txt"
tab=as.data.frame(dapc_noprior$ind.coord)
tab$INDIVIDUALS <- row.names(tab)
dpca_geo <- merge(x = tab, y=pop, by=c("INDIVIDUALS"))

### Add location information
SEBASTE_dpca$STRATA <- substr(SEBASTE_dpca$INDV, 1,5)

### Analyse how much percent of genetic variance is explained by each axe
pdf("Percent_dapc_all_loci.pdf")
percent= dapc2$eig/sum(dapc2$eig)*100
barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", names.arg=round(percent,2))
dev.off()
percent

#### Make a ggplot graph representing the PCA for the first and second axes considering the sampling locations
g1 = ggplot(SEBASTE_dpca, aes(x=DPC1, y=DPC2, fill=STRATA))+ geom_point(shape=21,size=3)+
  scale_fill_manual(values=rainbow(28))+
  labs(x="DPC1 (85.6%)")+
  labs(y="PC2 (6.0%)")+
  theme_classic()
g1


########### 4 - ACP  #################
results_acp =as.data.frame(dapc1$tab)
names(results_acp)
results_acp$INDV <- rownames(results_acp)
results_acp$STRATA <- substr(results_acp$INDV,1,5)
  
g1 = ggplot(results_acp, aes(x=`PCA-pc.1`, y=`PCA-pc.2`, fill=STRATA))+ geom_point(shape=21,size=3)+
  scale_fill_manual(values=rainbow(28))+
  labs(x="PC1")+
  labs(y="PC2")+
  theme_classic()
g1

### Trouver dans l'objet results_acp le pourcentage de variation expliqué par chaque axe.

########### 5 - DAPC avec prior #################

### Ajouter au fichier genind le nom des localités de chaque échantillon
SEBASTE_genind@pop <- pop$STRATA

### Vérifier que les localités ont bien été ajoutées au fichier genind
SEBASTE_genind@pop

### Find optimal alpha value
dapc_a_score <- dapc(SEBASTE_genind, n.da=100,n.pca=80)
temp_score <- optim.a.score(dapc_a_score)
# 13 Pcs to retain

## Run the DAPC with priors (prior as locations)
dapc_prior<-dapc(SEBASTE_genind,n.pca=13)
# 2 discriminant functions

### Visualize quickly the DPCA
myCol <- rainbow(28)
scatter(dapc_prior, posi.da="bottomright", bg="white",clabel = 0.3)

### Read the documentation for DPCA
?dapc

