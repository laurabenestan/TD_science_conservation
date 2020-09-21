# TD_science_conservation
TD pour l'UE Science de la conservation

---
title: "Script-dapc"
author: "Laura Benestan"
date: "9/19/2020"

---

# Analyse Discriminante en Composantes Principales (DACP)

Ce tutoriel fournit des recommandations pour appliquer l'analyse discriminante en composantes principales (DACP) sur un vaste ensemble de données génomiques telles que les données issues de RAD-sequencing. 

L'DACP agit comme une analyse en composantes principales (ACP) et accordera plus d'importance aux différences entre groupes qu'aux diff??rences intragroupes. Les groupes peuvent peuvent être spécifiés au début de l'analyse de manière aveugle (ci-dessous nommée DACP sans a priori) ou sur la base d'informations d'échantillonnage (DACP avec a priori). 

L'analyse DACP est particulièrement pertinente dans le cas des espèces marines, où l'on constate souvent une faible différenciation génétique qui empèche de révéler des différenecs génétiques subtiles mais significatives.

Cette approche statistique a été développée par Jombart et ses collègues par le biais du paquet "adegenet". Voir l'[excellent tutoriel](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf) réalisé par Thibault Jombart et Cattlin Collins.

## 0 - Préparez votre environnement R

- [R Version 3.6.0](https://cran.r-project.org/)
	
```{r}
library(radiator)
library(ggplot2)
library(adegenet)
library(ade4)
library(vcfR)
library(mapdata)
library(viridis)
```

## 1 - Télécharger des ensembles de données génomiques et spatiales

Importez le jeu de données en utilisant le paquet `vcfR`.
```{r}
vcf_sebastes <- read.vcfR("100snps_861ind_sebastes.recode.vcf")
SEBASTE_genind <- vcfR2genind(vcf_sebastes)
```

Vérifier le nombre d'échantillons.
```{r}
SEBASTE_genind
```

Importez le fichier contenant les informations suivantes : ID, LATITUDE, LONGITUDE et LIEU D'ÉCHANTILLONNAGE.
```{r}
pop <- read.table("population_map_sebastes.txt", header=TRUE, sep="\t", dec=".")
```

## 2 - Trouvez le nombre de clusters dans votre ensemble de données

Trouvez des clusters dans vos données en utilisant un objet genind.
```{r}
grp <- find.clusters(SEBASTE_genind, max.n.clust=10)
```

Choisissez le nombre d'axes (q) à conserver par rapport au nombre d'échantillons (N) contenus dans votre jeu de données et le nombre de variables utilisés (p). 
Ce nombre ne doit pas dépasser N/3 et q < p.
N'oubliez pas que plus vous sélectionnez d'axes (q), plus vous serez en mesure de détecter les regroupements dans les données.

Choisissez le nombre de groupes génétiques (K) optimal est celui qui correspond la plus petite valeur du critère d'information bayésien (BIC).
Le graphique vous montre le taux d'erreur minimum de validation croisée donné par le BIC. 
Ici, un K égal à 2-4 semble être le meilleur K à conserver. Prenez K = 4.

Observer la variation de valeur d'AIC en fonction du nombre de groupes.
```{r}
grp$Kstat
```

Vérifiez la taille du groupe pour chaque groupe génétique trouvé.
Cette étape vous donnera une idée de la validité de chaque groupe génétique trouvé.
```{r}
grp$size
```

Regarder quels individus appartiennent à quel groupe génétique et est ce que ça fait du sens?
```{r}
individuals <- head(grp$grp, 861) 
individuals
```

Sauvegarder l'information pour chaque individus
```{r}
data_kmeans <- data.frame(grp$grp)
data_kmeans$INDIVIDUALS <- row.names(data_kmeans) 
```

## 3 - Effectuer la DACP sans aucun a priori (sur la base du cluster K trouvé par l'analyse BIC)

Exécutez la DACP en utilisant au préalable le nombre d'axes égal à N/3 en considérant q = 1/2p

```{r}
dapc_noprior <-dapc(SEBASTE_genind, grp$grp)
```

Visualisez rapidement les résultats de la DACP.
```{r}
scatter(dapc_noprior, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
```

Explorez la contribution de chaque SNP à la discrimination des groupes génétiques trouvés.
```{r}
names(dapc_noprior)
head(dapc_noprior$var.contr)
load_dapc_noprior <- as.data.frame(dapc_noprior$var.contr)
```

Analyser le pourcentage de variation génétique observée à chacun des axes de l'analyse discriminante.
```{r}
pdf("Percent_dapc_all_loci.pdf")
percent= dapc_noprior$eig/sum(dapc_noprior$eig)*100
barplot(percent, ylab="Pourcentage de variation génétique expliquée par les eigenvectors", names.arg=round(percent,2))
dev.off()
percent
```

## 4 - Visualisez les résultats de la DACP avec la librairie ggplot

Sauvegardez les résultats de la DACP dans un tableau de données.
```{r}
tab=as.data.frame(dapc_noprior$ind.coord)
```

Ajoutez les informations sur chaque individu.
```{r}
tab$INDIVIDUALS <- row.names(tab)
```

Ajoutez l'information géographique.
```{r}
dpca_geo <- merge(x = tab, y=pop, by=c("INDIVIDUALS"))
```

Faire un graphique ggplot représentant les DACP pour le premier et le deuxième axe pour les localités.
Indiquez le pourcentage de variation expliqué par chaque axe.
```{r}
g = ggplot(dpca_geo, aes(x=LD1, y=LD2, fill=STRATA))+ geom_point(size=2, pch=21)+
  scale_fill_manual(values=rainbow(28))+
  labs(x="DPC1 (??)")+
  labs(y="DPC2 (??)")+
  theme_classic()
g
```

Sauvez le graphique ggplot.
```{r}
ggsave("DPCA_noprior.pdf",width=13,height=9,dpi=600,units="cm",useDingbats=F)
```

## 5 - Analysez des résultats d'ACP  #################

Sauvegardez les résultats de l'ACP à partir de l'information `tab`
```{r}
results_acp =as.data.frame(dapc_noprior$tab)
names(results_acp)
results_acp$INDV <- rownames(results_acp)
results_acp$STRATA <- substr(results_acp$INDV,1,5)
results_acp$SPECIES <- substr(results_acp$INDV,1,3)
```

Représentez les résultats avec la librarie ggplot
```{r}
g1 = ggplot(results_acp, aes(x=`PCA-pc.1`, y=`PCA-pc.2`, fill=STRATA))+ geom_point(shape=21,size=3)+
  scale_fill_manual(values=rainbow(28))+
  labs(x="PC1")+
  labs(y="PC2")+
  theme_classic()
g1
```

## 6 - Effectuez un DPCA avec a priori en se basant sur vos localités d'échantillonnage 

Trouver la valeur alpha optimale, qui représente le compromis entre le pouvoir de discrimination et le sur-ajustement. En cas de surajustement, vous pouvez observer un regroupement, même s'il s'agit de groupes aléatoires.

Ajoutez l'information des localités.
```{r}
SEBASTE_genind@pop <- pop$STRATA
```

Vérifier que les localités ont bien été ajoutées au fichier genind
```{r}
SEBASTE_genind@pop
```

Trouver la valeur alpha optimale c'est à dire le compromis entre le pouvoir de discrimination et le sur-ajustement.
Le alpha score est simplement la différence entre la proportion de réaffectation réussie par l'analyse (discrimination observée) et les valeurs obtenues en utilisant des groupes aléatoires (discrimination aléatoire). 
```{r}
dapc_a_score <- dapc(SEBASTE_genind, n.da=100,n.pca=80)
temp_score <- optim.a.score(dapc_a_score)
```

Exécutez la DACP avec le nombre d'axes indiqué par le alpha score. 
Si vous conservez trop d'axes, vous surchargez vos données, c'est pourquoi cette étape est vraiment cruciale.
```{r}
dapc_prior<-dapc(SEBASTE_genind,n.pca=13)
```

Visualisez rapidement les résultats de la DACP.
```{r}
scatter(dapc_prior, bg="white", scree.da=FALSE, legend=TRUE, solid=.4, col=rainbow(28))
```

Lisez la documentation pour améliorer votre graphique.
```{r}
?dapc
```

Faire un graphique ggplot de ces résutats en considérant les populations et les espèces.

## 6 - Represent the genetic clusters observed on a sampling map

Faire correspondre les coordonnées géographiques de chaque individu à son groupe génétique inféré afin de voir s'il y a un regroupement par rapport au lieu d'échantillonnage.
```{r}
kmean_geo <- merge(data_kmeans, pop, by="INDIVIDUALS")
```

Explorez les donnees geographiques.
```{r}
summary(pop)
```

Download the background data required for creating a map.
```{r}
wH <- map_data("worldHires", xlim=c(-72,-45), ylim=c(40,70)) # subset polygons surrounding 
```

```{r}
kmean_geo$Latitude <-as.numeric(as.character(kmean_geo$Latitude))
kmean_geo$Longitude<-as.numeric(as.character(kmean_geo$Longitude))
```

Make a map using `ggplot` and `sf` package.
```{r}
x_title="Longitude"
y_title="Latitude"
graph2 <- ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
  coord_fixed(xlim=c(-72,-45), ylim=c(40,70), ratio=1.2)+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  geom_point(aes(x = kmean_geo$Longitude, y = kmean_geo$Latitude, fill=kmean_geo$CLUSTER), shape=21)+
  labs(y=y_title)+  
  labs(x=x_title)+
  scale_fill_viridis(discrete=TRUE)+
theme_classic()
graph2
```

Save the graph.
```{r}
ggsave("Kmean_sebastes.pdf")
```
