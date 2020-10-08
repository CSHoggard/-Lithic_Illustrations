### EVALUATING THE SUITABILITY OF LITHIC ILLUSTRATIONS IN MORPHOMETRIC ANALYSES ###
### AUTHORS: CHRISTIAN STEVEN HOGGARD, THOMAS BIRCH, CORY MARIE STADE, KATRIEN JANIN AND FELIX RIEDE ###

### R SCRIPT AUTHOR: CHRISTIAN STEVEN HOGGARD ###
### EMAIL CONTACT: C.S.HOGGARD@SOTON.AC.UK ###
### OSF PAGE: https://osf.io/xtghn/ ###
### GITHUB: https://github.com/CSHoggard/-Lithic_Illustrations ###
### LAST EDITED: 14/04/2020 ###

### Abstract ###

### Illustrations of lithic artefacts are an abundant source of morphological and technological information for those interested in
### our human past. As a typical part of archaeological reports and publications, lithic drawings are - or have to be – trusted as 
### faithful reproductions of the selected artefacts. Despite the considerable epistemic work lithic illustrations (and illustrators)
### are expected to do, usually little information is available regarding the illustrator's technical skill; thus, it remains unknown
### whether drawings produced by illustrators of differing technical skill are comparable or produce images of equal analytical 
### potential to other media, e.g. photographs. The issue of lithic illustration accuracy is brought to the fore by the recent
### mergence of geometric morphometric approaches as innovative and powerful ways of describing and analysing complex shapes, as
### lithic illustrations provide one of the key sources for such analyses. Motivated by these issues, we present an experiment 
### investigating the degree of error observed in illustrations of differing technical illustrative skill. Analyses suggest that 
### lithic illustrations produced by individuals with a variety of experience in drawing lithics create, in the majority of instances,
### equally faithful representations (in outline shape) of chipped stone artefacts. With error observed in a small number of 
### instances, archaeologists are still urged to be critical of an illustration's source prior to lineal and geometric morphometric 
### methodologies. Despite this, archaeologists can be confident in their exactitude and we remain strong advocates in favour of 
### lithic illustrations as a readily available legacy resource for morphometric analyses. 

### r session information ###

### R version 3.4.3 (2017-11-30)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows >= 8 x64 (build 9200)

### base packages used: stats / graphics / grDevices utils / datasets / methods / base 
### attached other packages: MASS / vegan / cowplot / lattice / permute / forcats / stringr / dplyr / purrr / readr / tidyr
### tibble / ggplot2 / tidyverse / geomorph/ rgl / RRPP / psych

### For more extensive notes on this script, please refer to the R Markdown
### In improving the legibility of this script, all following additional information below is lower-case where possible

### stage 1: data download and package installation ###

setwd() ### set working directory
if(!require("psych")) install.packages('psych', repos='http://cran.us.r-project.org') ### psych 1.8.12
if(!require("geomorph")) install.packages('geomorph', repos='http://cran.us.r-project.org') ### geomorph 3.1.2
if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org') ### tidyverse 1.2.1
if(!require("vegan")) install.packages('vegan', repos='http://cran.us.r-project.org') ### vegan 2.5-4
if(!require("MASS")) install.packages('MASS', repos='http://cran.us.r-project.org') ### MASS 7.3-51.4
if(!require("cowplot")) install.packages('cowplot', repos='http://cran.us.r-project.org') ### cowplot v.0.9.3
if(!require("ggpubr")) install.packages('ggpubr', repos='http://cran.us.r-project.org') ### ggpubr v.0.2.4
if(!require("rio")) install.packages('ggpubr', repos='http://cran.us.r-project.org') ### rio v.0.5.16
if(!require("devtools")) install.packages('devtools', repos='http://cran.us.r-project.org') ### devtools v.2.3.1
devtools::install_github("akiopteryx/lambda") ### LaMBDA v.0.1.0.9000

library(psych) ### load the listed package
library(geomorph) ### load the listed package
library(tidyverse) ### load the listed package
library(vegan) ### load the listed package
library(MASS) ### load the listed package
library(cowplot) ### load the listed package
library(ggpubr) ### load the listed package
library(LaMBDA) ### load the listed package
library(rio) ### load the listed package

landmarks_elongated <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/landmarks_elongated.rds")
landmarks_handaxe <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/landmarks_handaxe.rds")
landmarks_tanged <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/landmarks_tanged.rds")
shape_data_elongated <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/shape_data_elongated.rds")
shape_data_tanged <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/shape_data_tanged.rds")
shape_data_handaxe <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/shape_data_handaxe.rds")
metric_data <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/metric_data.rds")
digitisation_error_landmarks <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/digitisation_error_landmarks.rds")
digitisation_error_landmarks_data <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/digitisation_error_landmarks_data.rds")
digitisation_error_metrics <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/digitisation_error_metrics.rds")
shape_data_sliders <- rio::import("https://github.com/CSHoggard/-Lithic_Illustrations/raw/master/shape_data_sliders.rds")

### stage 2: measuring intra-observor error ###

### stage 2a: digitisation and landmark error ###

gpa_digi_error <- gpagen(digitisation_error_landmarks, Proj = TRUE, curves =  shape_data_sliders, ProcD = TRUE, surfaces = NULL) ### generalised procrustes analysis: digitisation error
gpa_digi_error ### gpa details (including mean shape)
plot(gpa_digi_error) ### plots procrustes coordinates

gpa_digi_error_df <- geomorph.data.frame(gpa_digi_error, attempt = digitisation_error_landmarks_data$Attempt)
gpaprocD <- procD.lm(coords ~ attempt, data = gpa_digi_error_df) ### anova (shape vs. individual)
summary(gpaprocD) ### summary (p: 0.9915)
gpaprocD$aov.table$SS[1]/gpaprocD$aov.table$SS[3]*100 ### error as a percentage (8.6%)

### stage 2b: measurement error and landmark count examination ###

head(digitisation_error_metrics)
statsl  <- describe(digitisation_error_metrics$Length_mm) ### descriptive statistics
statsw  <- describe(digitisation_error_metrics$Width_mm) ### descriptive statistics
statssf <- describe(digitisation_error_metrics$Scale_Factor) ### descriptive statistics

(statsl$se/statsl$mean)   * 100 ### fractional uncertainty (length)
(statsw$se/statsw$mean)   * 100 ### fractional uncertainty (width)
(statssf$se/statssf$mean) * 100 ### fractional uncertainty (scale factor)

lasec(two.d.array(landmarks_elongated), 2, iter = 500) # may take some time
lasec(two.d.array(landmarks_tanged), 2, iter = 500) # may take some time
lasec(two.d.array(landmarks_handaxe), 2, iter = 500) # may take some time

### stage 3: generalised procrustes analysis (gpa) ###

gpa_elongated <- gpagen(landmarks_elongated, Proj = TRUE, ProcD = TRUE, curves = shape_data_sliders, surfaces = NULL) ### gpa: elongated examples (procrustes distance criterion)
gpa_elongated$coords ### procrustes coordinates (elongated examples)
plot(gpa_elongated) ### plot: procrustes coordinates (elongated examples)

gpa_tanged    <- gpagen(landmarks_tanged, Proj = TRUE, ProcD = TRUE, curves = shape_data_sliders, surfaces = NULL) ### gpa: tanged examples (procrustes distance criterion)
gpa_tanged$coords ### procrustes coordinates (tanged examples)
plot(gpa_tanged) ### plot: procrustes coordinates (tanged examples)

gpa_handaxe   <- gpagen(landmarks_handaxe, Proj = TRUE, ProcD = TRUE, curves = shape_data_sliders, surfaces = NULL) ### gpa: handaxe examples (procrustes distance criterion)
gpa_handaxe$coords ### procrustes coordinates (handaxe examples)
plot(gpa_handaxe) ### plot: procrustes coordinates (handaxe examples)

### stage 4: exploratory analysis of shape data (pca) ###

pca_elongated <- gm.prcomp(gpa_elongated$coords) ### pca (geomorph)
pca_elongated ### pca summary
elongated_ds <- cbind(shape_data_elongated, pca_elongated$x) ### tidyverse compatible format

pca_tanged <- gm.prcomp(gpa_tanged$coords) ### pca (geomorph)
pca_tanged ### pca summary
tanged_ds <- cbind(shape_data_tanged, pca_tanged$x) ### tidyverse compatible format

pca_handaxe <- gm.prcomp(gpa_handaxe$coords) ### pca (geomorph)
pca_handaxe ### pca summary
handaxe_ds <- cbind(shape_data_handaxe, pca_handaxe$x) ### tidyverse compatible format

figure_5 <- ggplot(data = elongated_ds) + geom_point(mapping = aes(x = PC1, y = PC2, colour = Artefact, shape = Class), size = 1) + labs(x = "Principal Component 1 (45.816%)", y = "Principal Component 2 (24.298%)", shape = "Method/skill", colour = "Artefact") + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")) +  scale_shape_manual(values=c(17,15,3,16)) + theme_minimal() + theme(text = element_text(size=8))
figure_6 <- ggplot(data = tanged_ds) + geom_point(mapping = aes(x = PC1, y = PC2, colour = Artefact, shape = Class), size = 1) + labs(x = "Principal Component 1 (50.805%)", y = "Principal Component 2 (33.926%)", shape = "Method/skill", colour = "Artefact") + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")) +  scale_shape_manual(values=c(17,15,3,16)) + theme_minimal() + theme(text = element_text(size=8))
figure_7 <- ggplot(data = handaxe_ds) + geom_point(mapping = aes(x = PC1, y = PC2, colour = Artefact, shape = Class), size = 1) + labs(x = "Principal Component 1 (84.840%)", y = "Principal Component 2 (6.987%)", shape = "Method/skill", colour = "Artefact") + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")) +  scale_shape_manual(values=c(17,15,3,16)) + theme_minimal() + theme(text = element_text(size=8))

ggsave("Fig5.tiff", plot = figure_5, dpi = 300, units = "mm", height = 90, width = 120)
ggsave("Fig6.tiff", plot = figure_6, dpi = 400, units = "mm", height = 90, width = 120)
ggsave("Fig7.tiff", plot = figure_7, dpi = 400, units = "mm", height = 90, width = 120)

### stage 5: analytical framework for the shape variables ###

### stage 5a: elongated examples ###

elongated1 <- gpa_elongated$coords[, , 1:18] 
df_elongated1 <- geomorph.data.frame(shape = elongated1, class = shape_data_elongated$Class[1:18], artefact = shape_data_elongated$Artefact[1:18]) 
E1 <- procD.lm(shape ~ class, data = df_elongated1, print.progress = FALSE) 
summary(E1) 

elongated1pcs <- as.data.frame(pca_elongated$x[1:18, 1:10])
elongated1class <- shape_data_elongated$Class[1:18] 
elongated1pcs <- cbind(elongated1pcs, elongated1class)
elongated1pcs <- rename(elongated1pcs, Class = elongated1class)
elongated1lda <- lda(Class ~ ., data = elongated1pcs)
elongated1ldapredict <- predict(elongated1lda) 
elongated1ldaplot <- cbind(elongated1ldapredict$x[, 1:3], elongated1pcs) 
elongated1ldaplotggplot <- ggplot(elongated1ldaplot, aes(LD1, LD2)) + 
  geom_point(aes(shape = Class), size = 2) + 
  labs(x = "LDA 1", y = "LDA 2") + 
  scale_shape_manual(values=c(17,15,3,16)) + 
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6), 
        legend.position = "none")

elongated2 <- gpa_elongated$coords[, , 19:36]
df_elongated2 <- geomorph.data.frame(shape = elongated2, class = shape_data_elongated$Class[19:36], artefact = shape_data_elongated$Artefact[19:36])
E2 <- procD.lm(shape ~ class, data = df_elongated2, print.progress = FALSE)
summary(E2)

elongated2pcs <- as.data.frame(pca_elongated$x[19:36, 1:10])
elongated2class <- shape_data_elongated$Class[19:36]
elongated2pcs <- cbind(elongated2pcs, elongated2class)
elongated2pcs <- rename(elongated2pcs, Class = elongated2class)
elongated2lda <- lda(Class ~ ., elongated2pcs)
elongated2ldapredict <- predict(elongated2lda)
elongated2ldaplot <- cbind(elongated2ldapredict$x[, 1:3], elongated2pcs)
elongated2ldaplotggplot <- ggplot(elongated2ldaplot, aes(LD1, LD2)) + 
  geom_point(aes(shape = Class), size = 2) + 
  labs(x = "LDA 1", y = "LDA 2") + 
  scale_shape_manual(values=c(17,15,3,16)) + 
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 6),
        axis.title.x = element_text(size = 6))

elongated3 <- gpa_elongated$coords[, , 37:54]
df_elongated3 <- geomorph.data.frame(shape = elongated3, class = shape_data_elongated$Class[37:54], artefact = shape_data_elongated$Artefact[37:54])
E3 <- procD.lm(shape ~ class, data = df_elongated3, print.progress = FALSE)
summary(E3) 

elongated3pcs <- as.data.frame(pca_elongated$x[37:54, 1:10])
elongated3class <- shape_data_elongated$Class[37:54]
elongated3pcs <- cbind(elongated3pcs, elongated3class)
elongated3pcs <- rename(elongated3pcs, Class = elongated3class)
elongated3lda <- lda(Class ~ ., elongated3pcs)
elongated3ldapredict <- predict(elongated3lda)
elongated3ldaplot <- cbind(elongated3ldapredict$x[, 1:3], elongated3pcs)
elongated3ldaplotggplot <- ggplot(elongated3ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

elongated4 <- gpa_elongated$coords[, , 55:72]
df_elongated4 <- geomorph.data.frame(shape = elongated4, class = shape_data_elongated$Class[55:72], artefact = shape_data_elongated$Artefact[55:72])
E4 <- procD.lm(shape ~ class, data = df_elongated4, print.progress = FALSE)
summary(E4)

elongated4pcs <- as.data.frame(pca_elongated$x[55:72, 1:10])
elongated4class <- shape_data_elongated$Class[55:72]
elongated4pcs <- cbind(elongated4pcs, elongated4class)
elongated4pcs <- rename(elongated4pcs, Class = elongated4class)
elongated4lda <- lda(Class ~ ., elongated4pcs)
elongated4ldapredict <- predict(elongated4lda)
elongated4ldaplot <- cbind(elongated4ldapredict$x[, 1:3], elongated4pcs)
elongated4ldaplotggplot <- ggplot(elongated4ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

elongated5 <- gpa_elongated$coords[, , 73:90] 
df_elongated5 <- geomorph.data.frame(shape = elongated5, class = shape_data_elongated$Class[73:90], artefact = shape_data_elongated$Artefact[73:90]) 
E5 <- procD.lm(shape ~ class, data = df_elongated5, print.progress = FALSE)
summary(E5)

elongated5pcs <- as.data.frame(pca_elongated$x[73:90, 1:10])
elongated5class <- shape_data_elongated$Class[73:90]
elongated5pcs <- cbind(elongated5pcs, elongated5class)
elongated5pcs <- rename(elongated5pcs, Class = elongated5class)
elongated5lda <- lda(Class ~ ., elongated5pcs)
elongated5ldapredict <- predict(elongated5lda)
elongated5ldaplot <- cbind(elongated5ldapredict$x[, 1:3], elongated5pcs)
elongated5ldaplotggplot <- ggplot(elongated5ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 6),
        axis.title.x = element_text(size = 6))

tanged1 <- gpa_tanged$coords[, , 1:18]
df_tanged1 <- geomorph.data.frame(shape = tanged1, class = shape_data_tanged$Class[1:18], artefact = shape_data_tanged$Artefact[1:18])
T1 <- procD.lm(shape ~ class, data = df_tanged1, print.progress = FALSE)
summary(T1) 

tanged1pcs <- as.data.frame(pca_tanged$x[1:18, 1:10])
tanged1class <- shape_data_tanged$Class[1:18]
tanged1pcs <- cbind(tanged1pcs, tanged1class)
tanged1pcs <- rename(tanged1pcs, Class = tanged1class)
tanged1lda <- lda(Class ~ ., tanged1pcs)
tanged1ldapredict <- predict(tanged1lda)
tanged1ldaplot <- cbind(tanged1ldapredict$x[, 1:3], tanged1pcs)
tanged1ldaplotggplot <- ggplot(tanged1ldaplot, aes(LD1, LD2)) + 
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

tanged2 <- gpa_tanged$coords[, , 19:36]
df_tanged2 <- geomorph.data.frame(shape = tanged2, class = shape_data_tanged$Class[19:36], artefact = shape_data_tanged$Artefact[19:36])
T2 <- procD.lm(shape ~ class, data = df_tanged2, print.progress = FALSE)
summary(T2)

tanged2pcs <- as.data.frame(pca_tanged$x[19:36, 1:10])
tanged2class <- shape_data_tanged$Class[19:36]
tanged2pcs <- cbind(tanged2pcs, tanged2class)
tanged2pcs <- rename(tanged2pcs, Class = tanged2class)
tanged2lda <- lda(Class ~ ., tanged2pcs)
tanged2ldapredict <- predict(tanged2lda)
tanged2ldaplot <- cbind(tanged2ldapredict$x[, 1:3], tanged2pcs)
tanged2ldaplotggplot <- ggplot(tanged2ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

tanged3 <- gpa_tanged$coords[, , 37:54] 
df_tanged3 <- geomorph.data.frame(shape = tanged3, class = shape_data_tanged$Class[37:54], artefact = shape_data_tanged$Artefact[37:54])
T3 <- procD.lm(shape ~ class, data = df_tanged3, print.progress = FALSE)
summary(T3)

tanged3pcs <- as.data.frame(pca_tanged$x[37:54, 1:10])
tanged3class <- shape_data_tanged$Class[37:54]
tanged3pcs <- cbind(tanged3pcs, tanged3class)
tanged3pcs <- rename(tanged3pcs, Class = tanged3class)
tanged3lda <- lda(Class ~ ., tanged3pcs)
tanged3ldapredict <- predict(tanged3lda)
tanged3ldaplot <- cbind(tanged3ldapredict$x[, 1:3], tanged3pcs)
tanged3ldaplotggplot <- ggplot(tanged3ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

tanged4 <- gpa_tanged$coords[, , 55:72]
df_tanged4 <- geomorph.data.frame(shape = tanged4, class = shape_data_tanged$Class[55:72], artefact = shape_data_tanged$Artefact[55:72])
T4 <- procD.lm(shape ~ class, data = df_tanged4, print.progress = FALSE)
summary(T4)

tanged4pcs <- as.data.frame(pca_tanged$x[55:72, 1:10])
tanged4class <- shape_data_tanged$Class[55:72]
tanged4pcs <- cbind(tanged4pcs, tanged4class)
tanged4pcs <- rename(tanged4pcs, Class = tanged4class)
tanged4lda <- lda(Class ~ ., tanged4pcs)
tanged4ldapredict <- predict(tanged4lda)
tanged4ldaplot <- cbind(tanged4ldapredict$x[, 1:3], tanged4pcs)
tanged4ldaplotggplot <- ggplot(tanged4ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

tanged5 <- gpa_tanged$coords[, , 73:90]
df_tanged5 <- geomorph.data.frame(shape = tanged5, class = shape_data_tanged$Class[73:90], artefact = shape_data_tanged$Artefact[73:90])
T5 <- procD.lm(shape ~ class, data = df_tanged5, print.progress = FALSE)
summary(T5)

tanged5pcs <- as.data.frame(pca_tanged$x[73:90, 1:10])
tanged5class <- shape_data_tanged$Class[73:90]
tanged5pcs <- cbind(tanged5pcs, tanged5class)
tanged5pcs <- rename(tanged5pcs, Class = tanged5class)
tanged5lda <- lda(Class ~ ., tanged5pcs)
tanged5ldapredict <- predict(tanged5lda)
tanged5ldaplot <- cbind(tanged5ldapredict$x[, 1:3], tanged5pcs)
tanged5ldaplotggplot <- ggplot(tanged5ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") + 
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

handaxe1 <- gpa_handaxe$coords[, , 1:18]
df_handaxe1 <- geomorph.data.frame(shape = handaxe1, class = shape_data_handaxe$Class[1:18], artefact = shape_data_handaxe$Artefact[1:18])
H1 <- procD.lm(shape ~ class, data = df_handaxe1, print.progress = FALSE)
summary(H1)

handaxe1pcs <- as.data.frame(pca_handaxe$x[1:18, 1:10])
handaxe1class <- shape_data_handaxe$Class[1:18]
handaxe1pcs <- cbind(handaxe1pcs, handaxe1class)
handaxe1pcs <- rename(handaxe1pcs, Class = handaxe1class)
handaxe1lda <- lda(Class ~ ., handaxe1pcs)
handaxe1ldapredict <- predict(handaxe1lda)
handaxe1ldaplot <- cbind(handaxe1ldapredict$x[, 1:3], handaxe1pcs)
handaxe1ldaplotggplot <- ggplot(handaxe1ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") + 
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

handaxe2 <- gpa_handaxe$coords[, , 19:36]
df_handaxe2 <- geomorph.data.frame(shape = handaxe2, class = shape_data_handaxe$Class[19:36], artefact = shape_data_handaxe$Artefact[19:36])
H2 <- procD.lm(shape ~ class, data = df_handaxe2, print.progress = FALSE)
summary(H2)

handaxe2pcs <- as.data.frame(pca_handaxe$x[19:36, 1:10]) 
handaxe2class <- shape_data_handaxe$Class[19:36]
handaxe2pcs <- cbind(handaxe2pcs, handaxe2class)
handaxe2pcs <- rename(handaxe2pcs, Class = handaxe2class)
handaxe2lda <- lda(Class ~ ., handaxe2pcs)
handaxe2ldapredict <- predict(handaxe2lda)
handaxe2ldaplot <- cbind(handaxe2ldapredict$x[, 1:3], handaxe2pcs)
handaxe2ldaplotggplot <- ggplot(handaxe2ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

handaxe3 <- gpa_handaxe$coords[, , 37:54]
df_handaxe3 <- geomorph.data.frame(shape = handaxe3, class = shape_data_handaxe$Class[37:54], artefact = shape_data_handaxe$Artefact[37:54])
H3 <- procD.lm(shape ~ class, data = df_handaxe3, print.progress = FALSE)
summary(H3)

handaxe3pcs <- as.data.frame(pca_handaxe$x[37:54, 1:10])
handaxe3class <- shape_data_handaxe$Class[37:54]
handaxe3pcs <- cbind(handaxe3pcs, handaxe3class)
handaxe3pcs <- rename(handaxe3pcs, Class = handaxe3class)
handaxe3lda <- lda(Class ~ ., handaxe3pcs)
handaxe3ldapredict <- predict(handaxe3lda)
handaxe3ldaplot <- cbind(handaxe3ldapredict$x[, 1:3], handaxe3pcs)
handaxe3ldaplotggplot <- ggplot(handaxe3ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") + 
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

handaxe4 <- gpa_handaxe$coords[, , 55:72] 
df_handaxe4 <- geomorph.data.frame(shape = handaxe4, class = shape_data_handaxe$Class[55:72], artefact = shape_data_handaxe$Artefact[55:72]) 
H4 <- procD.lm(shape ~ class, data = df_handaxe4, print.progress = FALSE) 
summary(H4)

handaxe4pcs <- as.data.frame(pca_handaxe$x[55:72, 1:10])
handaxe4class <- shape_data_handaxe$Class[55:72]
handaxe4pcs <- cbind(handaxe4pcs, handaxe4class)
handaxe4pcs <- rename(handaxe4pcs, Class = handaxe4class)
handaxe4lda <- lda(Class ~ ., handaxe4pcs)
handaxe4ldapredict <- predict(handaxe4lda)
handaxe4ldaplot <- cbind(handaxe4ldapredict$x[, 1:3], handaxe4pcs)
handaxe4ldaplotggplot <- ggplot(handaxe4ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

handaxe5 <- gpa_handaxe$coords[, , 73:90]
df_handaxe5 <- geomorph.data.frame(shape = handaxe5, class = shape_data_handaxe$Class[73:90], artefact = shape_data_handaxe$Artefact[73:90])
H5 <- procD.lm(shape ~ class, data = df_handaxe5, print.progress = FALSE)
summary(H5)

handaxe5pcs <- as.data.frame(pca_handaxe$x[73:90, 1:10])
handaxe5class <- shape_data_handaxe$Class[73:90]
handaxe5pcs <- cbind(handaxe5pcs, handaxe5class)
handaxe5pcs <- rename(handaxe5pcs, Class = handaxe5class)
handaxe5lda <- lda(Class ~ ., handaxe5pcs)
handaxe5ldapredict <- predict(handaxe5lda)
handaxe5ldaplot <- cbind(handaxe5ldapredict$x[, 1:3], handaxe5pcs)
handaxe5ldaplotggplot <- ggplot(handaxe5ldaplot, aes(LD1, LD2)) +
  geom_point(aes(shape = Class), size = 2) +
  labs(x = "LDA 1", y = "LDA 2") +
  scale_shape_manual(values=c(17,15,3,16)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none",  
        axis.title.y = element_text(size = 6),  
        axis.title.x = element_text(size = 6))

lda.figure <-  plot_grid(elongated1ldaplotggplot,
                         elongated2ldaplotggplot,
                         elongated3ldaplotggplot, 
                         elongated4ldaplotggplot, 
                         elongated5ldaplotggplot,
                         tanged1ldaplotggplot, 
                         tanged2ldaplotggplot, 
                         tanged3ldaplotggplot, 
                         tanged4ldaplotggplot, 
                         tanged5ldaplotggplot,
                         handaxe1ldaplotggplot,
                         handaxe2ldaplotggplot,
                         handaxe3ldaplotggplot,
                         handaxe4ldaplotggplot,
                         handaxe5ldaplotggplot,
                         ncol = 5, align = 'v')

lda.figure.legend <- get_legend(elongated1ldaplotggplot +
                                  guides(color = guide_legend(nrow = 1)) +
                                  theme(legend.position = "bottom"))

plot_grid(lda.figure, lda.figure.legend, ncol = 1, rel_heights = c(1, .1))

### stage 6: exploratory framework for the measurement data ###

metric_data_elongated <- metric_data[which(metric_data$Type=="Elongated"),]
metric_data_handaxe   <- metric_data[which(metric_data$Type=="Handaxe"),]
metric_data_tanged    <- metric_data[which(metric_data$Type=="Tanged"),]

figure_8a <- ggplot(metric_data_elongated, aes(Length_mm, Width_mm, colour = Artefact, shape = Class)) +
  geom_point() + 
  facet_grid(cols = vars(Artefact), scales = "free", labeller=label_both) + 
  labs(x = "Length (mm)", y = "Width (mm)") + 
  scale_shape_manual(values=c(17,15,3,16)) + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none")

figure_8b <- ggplot(metric_data_tanged, aes(Length_mm, Width_mm, colour = Artefact, shape = Class)) + 
  geom_point() + 
  facet_grid(cols = vars(Artefact), scales = "free", labeller=label_both) + 
  labs(x = "Length (mm)", y = "Width (mm)") +  
  scale_shape_manual(values=c(17,15,3,16)) + 
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none")

figure_8c <- ggplot(metric_data_handaxe, aes(Length_mm, Width_mm, colour = Artefact, shape = Class)) +
  geom_point() + 
  facet_grid(cols = vars(Artefact), scales = "free", labeller=label_both) + 
  labs(x = "Length (mm)", y = "Width (mm)") +  
  scale_shape_manual(values=c(17,15,3,16)) +  
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none")

figure_8 <- plot_grid(figure_8a, 
                      figure_8b, 
                      figure_8c, 
                      labels= "AUTO", 
                      ncol = 1, 
                      align = 'v') 

plot_grid(figure_8, lda.figure.legend, ncol = 1, rel_heights = c(1, .1))

ggsave("Fig8.tiff", plot = last_plot(), dpi = 300, units = "mm", height = 150, width = 130)

### stage 7: analytical framework for the measurement data ###

### stage 7a: elongated examples ###

metric_data_elongated <- arrange(metric_data_elongated, Artefact, Class)

metric_data_elongated[1:18,]  %>% 
  group_by(Class) %>% 
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[1:18,]  %>% 
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% 
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[19:36,] %>% 
  group_by(Class) %>% 
  summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[19:36,]  %>% 
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[37:54,] %>% 
  group_by(Class) %>% 
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[37:54,]  %>% 
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% 
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[55:72,] %>% 
  group_by(Class) %>% 
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[55:72,]  %>% 
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% 
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[73:90,] %>% 
  group_by(Class) %>% 
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[73:90,]  %>% 
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% 
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

shapiro.test(metric_data_elongated$Length_mm[2:4])
shapiro.test(metric_data_elongated$Width_mm[2:4])
shapiro.test(metric_data_elongated$Length_mm[5:17]) 
shapiro.test(metric_data_elongated$Width_mm[5:17])
shapiro.test(metric_data_elongated$Length_mm[20:22])
shapiro.test(metric_data_elongated$Width_mm[20:22]) 
shapiro.test(metric_data_elongated$Length_mm[23:35])
shapiro.test(metric_data_elongated$Width_mm[23:35])
shapiro.test(metric_data_elongated$Length_mm[38:40]) 
shapiro.test(metric_data_elongated$Width_mm[38:40]) 
shapiro.test(metric_data_elongated$Length_mm[41:53])
shapiro.test(metric_data_elongated$Width_mm[41:53])
shapiro.test(metric_data_elongated$Length_mm[56:58])
shapiro.test(metric_data_elongated$Width_mm[56:58])
shapiro.test(metric_data_elongated$Length_mm[59:71])
shapiro.test(metric_data_elongated$Width_mm[59:71])
shapiro.test(metric_data_elongated$Length_mm[74:76])
shapiro.test(metric_data_elongated$Width_mm[74:76])
shapiro.test(metric_data_elongated$Length_mm[77:89])
shapiro.test(metric_data_elongated$Width_mm[77:89])

summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_elongated[1:18,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_elongated[1:18,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_elongated[19:36,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_elongated[19:36,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_elongated[37:54,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_elongated[37:54,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_elongated[55:72,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_elongated[55:72,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_elongated[73:90,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_elongated[73:90,]))

metric_data_handaxe <- arrange(metric_data_handaxe, Artefact, Class)

metric_data_handaxe[1:18,]  %>% 
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% 
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[1:18,]  %>% 
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% 
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[19:36,] %>%
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[19:36,]  %>% 
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[37:54,] %>% 
  group_by(Class) %>%
  summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[37:54,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[55:72,] %>%
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[55:72,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[73:90,] %>%
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE),
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[73:90,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

shapiro.test(metric_data_handaxe$Length_mm[2:4])
shapiro.test(metric_data_handaxe$Width_mm[2:4])
shapiro.test(metric_data_handaxe$Length_mm[5:17])
shapiro.test(metric_data_handaxe$Width_mm[5:17])
shapiro.test(metric_data_handaxe$Length_mm[20:22])
shapiro.test(metric_data_handaxe$Width_mm[20:22])
shapiro.test(metric_data_handaxe$Length_mm[23:35])
shapiro.test(metric_data_handaxe$Width_mm[23:35])
shapiro.test(metric_data_handaxe$Length_mm[38:40])
shapiro.test(metric_data_handaxe$Width_mm[38:40])
shapiro.test(metric_data_handaxe$Length_mm[41:53])
shapiro.test(metric_data_handaxe$Width_mm[41:53])
shapiro.test(metric_data_handaxe$Length_mm[56:58])
shapiro.test(metric_data_handaxe$Width_mm[56:58])
shapiro.test(metric_data_handaxe$Length_mm[59:71])
shapiro.test(metric_data_handaxe$Width_mm[59:71]) 
shapiro.test(metric_data_handaxe$Length_mm[74:76])
shapiro.test(metric_data_handaxe$Width_mm[74:76])
shapiro.test(metric_data_handaxe$Length_mm[77:89])
shapiro.test(metric_data_handaxe$Width_mm[77:89])

summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_handaxe[1:18,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_handaxe[1:18,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_handaxe[19:36,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_handaxe[19:36,]))
adonis(cbind(Length_mm, Width_mm) ~ Class, 
       data = metric_data_handaxe[37:54,])
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_handaxe[55:72,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_handaxe[55:72,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_handaxe[73:90,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_handaxe[73:90,]))

metric_data_tanged <- arrange(metric_data_tanged, Artefact, Class)

metric_data_tanged[1:18,]  %>% 
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[1:18,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[19:36,] %>%
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[19:36,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[37:54,] %>%
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[37:54,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[55:72,] %>%
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[55:72,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[73:90,] %>%
  group_by(Class) %>%
  summarise(count = n(), 
            mean = mean(Length_mm, na.rm = TRUE), 
            sd = sd(Length_mm, na.rm = TRUE), 
            min = min(Length_mm, na.rm = TRUE), 
            max = max(Length_mm, na.rm = TRUE), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[73:90,]  %>%
  filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>%
  summarise(count = n(), 
            cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>%
  dplyr::mutate_if(is.numeric, format, 1)

shapiro.test(metric_data_tanged$Length_mm[2:4]) 
shapiro.test(metric_data_tanged$Width_mm[2:4])
shapiro.test(metric_data_tanged$Length_mm[5:17])
shapiro.test(metric_data_tanged$Width_mm[5:17])
shapiro.test(metric_data_tanged$Length_mm[20:22])
shapiro.test(metric_data_tanged$Width_mm[20:22])
shapiro.test(metric_data_tanged$Length_mm[23:35])
shapiro.test(metric_data_tanged$Width_mm[23:35])
shapiro.test(metric_data_tanged$Length_mm[38:40])
shapiro.test(metric_data_tanged$Width_mm[38:40])
shapiro.test(metric_data_tanged$Length_mm[41:53])
shapiro.test(metric_data_tanged$Width_mm[41:53])
shapiro.test(metric_data_tanged$Length_mm[56:58])
shapiro.test(metric_data_tanged$Width_mm[56:58])
shapiro.test(metric_data_tanged$Length_mm[59:71])
shapiro.test(metric_data_tanged$Width_mm[59:71])
shapiro.test(metric_data_tanged$Length_mm[74:76]) 
shapiro.test(metric_data_tanged$Width_mm[74:76]) 
shapiro.test(metric_data_tanged$Length_mm[77:89]) 
shapiro.test(metric_data_tanged$Width_mm[77:89]) 

summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_tanged[1:18,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_tanged[1:18,])) 
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_tanged[19:36,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_tanged[19:36,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_tanged[37:54,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_tanged[37:54,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_tanged[55:72,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_tanged[55:72,]))
summary(manova(cbind(Length_mm, Width_mm) ~ Class, 
               data = metric_data_tanged[73:90,]))
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, 
                   data = metric_data_tanged[73:90,]))

########
