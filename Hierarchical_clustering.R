library(vegan)
library(cluster)
library(beepr)
library(dplyr)
library(ggplot2)

data("varechem")
chem_d <- select(varechem, N:Mo)
chem_d
varechem %>%
  select(N:Mo) %>%
  pairs(upper.panel = NULL)

chem_m <- round(vegdist(chem_d, method = "euclidean"), 1)

length(chem_m)
as.matrix(chem_m)[c(1:5), c(1:5)]

plot(chem_clust, las = 1, 
     main="Cluster diagram of soil chemistry", 
     xlab="Sample", ylab="Euclidean distance")
rect.hclust(chem_clust, 2, border="red") 
rect.hclust(chem_clust, 4, border="blue")
rect.hclust(chem_clust, 5, border="darkgreen")

fit <- cascadeKM(chem_d, 1, 10, iter = 5000)

plot(fit, sortg = TRUE, grpmts.plot = TRUE)

fit$results %>% 
  as.data.frame() %>%
  rownames_to_column("metric") %>%
  pivot_longer(names_to = "groups", 
               values_to = "value", 
               - metric) %>%
  mutate(groups = str_extract(groups, "\\d+"), 
         groups = as.numeric(groups)) %>%
  filter(metric != "SSE") %>%
  ggplot(aes(x=groups, y = value)) + theme_bw(16) +
  geom_line(lwd=1.5, col="blue") +
  geom_point(pch=21, col="lightgrey", 
             bg="blue", stroke = 1.5, size=5) +
  scale_x_continuous(breaks = c(2:10), labels = c(2:10)) +
  theme(panel.grid.minor.x = element_blank()) 

grps <- as_tibble(fit$partition)
grps

plot(chem_clust, las = 1,
     label = grps$`5 groups`,
     main = "Cluster diagram of soil chemistry",
     xlab = "Sample",
     ylab = "Euclidean distance")
rect.hclust(chem_clust, 5, border = "red")



#Try using Caatinga data
df <- read.csv("L:/N/Lucy_Caatinga/overall_fuzzy_df.csv");beep()
names(df)

df2 <- df[,c('x', 'y', 'PCA1', 'PCA2', 'PCA3')] #select PCA columns, and also the x and y columns so that you know which pixels are selected
head(df2)
df2 <- na.omit(df2)

#df2 %>%
 # pairs(upper.panel = NULL) #takes a long time have saved plot

a <- dplyr::sample_frac(df2, 0.01)

#df2_m <- vegdist(dplyr::sample_frac(df2, method = "euclidean")
df2_m <- vegdist(a, method = "euclidean")

hc <- hclust(df2_m)

plot(hc)
rect.hclust(hc , k = 8, border = 2:6)
abline(h = 8, col = 'red')

clusters = 6
groups <- cutree(hc, h = 8)
groups
groups <- as.data.frame(groups)
groups <- merge(a, groups)

clusplot(a, groups, color=TRUE, shade=TRUE, labels=4, lines=0);beep()

plot(hc)

#try bootstrapping to test cluster adequacy
library(fpc)
clus.boot <- clusterboot(a, 
                         B=100, # Number of bootstrap resamples
                         clustermethod=hclustCBI, # for hierarchical clustering 
                         method="ward.D", # use what we used in "hclust"
                         k=clusters, 
                         count=FALSE);beep() # Show progress on screen?
clus.boot
AvgJaccard <- clus.boot$bootmean
Instability <- clus.boot$bootbrd/1000
Clusters <- c(1:clusters)
Eval <- cbind(Clusters, AvgJaccard, Instability)
Eval


#how about on the raw veg param data?
names(df)

dfveg <- df[,c('biomass_ne', 'height_ne', 'LAI_dry_season_ne', 'LAI_wet_season_ne', 'sea_ne', 'biomass_err_ne', 'fire_log')] #select PCA columns 
head(dfveg)
dfveg <- na.omit(dfveg)

a_veg <- dplyr::sample_frac(dfveg, 0.01)

dfveg_m <- vegdist(a_veg, method = "euclidean")
hc_veg <- hclust(dfveg_m)

plot(hc_veg)
clusplot(a_veg, groups, color=TRUE, shade=TRUE, labels=4, lines=0);beep()

#need to sample with the x and y coordinates, in order to know what the env variables are at the sampled pixels


#multiple sampling https://infer.netlify.app/reference/rep_sample_n.html 

t <- rep_slice_sample(df2,
                      prop = 0.01,
                      replace = FALSE, 
                      weight_by = NULL,
                      reps = 10)
slices <- df2 %>%
  rep_slice_sample(prop = 0.01, reps = 10)

library(infer)

