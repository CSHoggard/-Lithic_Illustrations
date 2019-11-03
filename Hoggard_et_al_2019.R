### EVALUATING THE SUITABILITY OF LITHIC ILLUSTRATIONS IN MORPHOMETRIC ANALYSES ###
### AUTHORS: CHRISTIAN STEVEN HOGGARD, THOMAS BIRCH, CORY MARIE STADE, KATRIEN JANIN AND FELIX RIEDE ###

### R SCRIPT AUTHOR: CHRISTIAN STEVEN HOGGARD ###
### EMAIL CONTACT: C.S.HOGGARD@SOTON.AC.UK ###
### OSF PAGE: https://osf.io/xtghn/ ###
### GITHUB: https://github.com/CSHoggard/-Lithic_Illustrations_AASC ###
### LAST EDITED: 01/11/2019 ###

### ABSTRACT ###

### Illustrations of lithic artefacts are an abundant source of morphological and technological information for those interested in
### our human past. As a typical part of archaeological reports and publications, lithic drawings are - or have to be – trusted as 
### faithful reproductions of the selected artefacts. Despite the considerable epistemic work lithic illustrations (and illustrators)
### are expected to do, usually little information is available regarding the illustrator’s technical skill; thus, it remains unknown
### whether drawings produced by illustrators of differing technical skill are comparable or produce images of equal analytical 
### potential to other media, e.g. photographs. The issue of lithic illustration accuracy is brought to the fore by the recent
### mergence of geometric morphometric approaches as innovative and powerful ways of describing and analysing complex shapes, as
### lithic illustrations provide one of the key sources for such analyses. Motivated by these issues, we present an experiment 
### investigating the degree of error observed in illustrations of differing technical illustrative skill. Analyses suggest that 
### lithic illustrations produced by individuals with a variety of experience in drawing lithics create, in the majority of instances,
### equally faithful representations (in outline shape) of chipped stone artefacts. With error observed in a small number of 
### instances, archaeologists are still urged to be critical of an illustration’s source prior to lineal and geometric morphometric 
### methodologies. Despite this, archaeologists can be confident in their exactitude and we remain strong advocates in favour of 
### lithic illustrations as a readily available legacy resource for morphometric analyses. 

### R SESSION INFORMATION ###

# R version 3.4.3 (2017-11-30)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)

# Base packages used: stats / graphics / grDevices utils / datasets / methods / base 
# Attached other packages: MASS / vegan / cowplot / lattice / permute / forcats / stringr / dplyr / purrr / readr / tidyr
# tibble / ggplot2 / tidyverse / geomorph/ rgl / RRPP / psych

# For more extensive notes on this script, please refer to the R Markdown

### STAGE 1: DATA COLLECTION AND PACKAGE INSTALLATION ###

setwd() # SET WORKING DIRECTORY
if(!require("psych")) install.packages('psych', repos='http://cran.us.r-project.org') #PSYCH 1.8.12
if(!require("geomorph")) install.packages('geomorph', repos='http://cran.us.r-project.org') #GEOMORPH 3.1.2
if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org') #TIDYVERSE 1.2.1
if(!require("vegan")) install.packages('vegan', repos='http://cran.us.r-project.org') #vegan 2.5-4
if(!require("MASS")) install.packages('MASS', repos='http://cran.us.r-project.org') #MASS 7.3-51.4
if(!require("cowplot")) install.packages('cowplot', repos='http://cran.us.r-project.org') #cowplot v.0.9.3
if(!require("ggpubr")) install.packages('cowplot', repos='http://cran.us.r-project.org') #cowplot v.0.9.3

library(psych) #LOAD THE LISTED PACKAGE
library(geomorph) #LOAD THE LISTED PACKAGE
library(tidyverse) #LOAD THE LISTED PACKAGE
library(vegan) #LOAD THE LISTED PACKAGE
library(MASS) #LOAD THE LISTED PACKAGE
library(cowplot) #LOAD THE LISTED PACKAGE
library(ggpubr) #LOAD THE LISTED PACKAGE

landmarks_elongated   <- readland.tps("elongated.tps", readcurves = TRUE) #LANDMARK DATA FOR ALL ELONGATED EXAMPLES (SAFE TO IGNORE RESCALING ISSUE)
landmarks_tanged      <- readland.tps("tanged.tps", readcurves = TRUE)    #LANDMARK DATA FOR ALL TANGED EXAMPLES (SAFE TO IGNORE RESCALING ISSUE)
landmarks_handaxe     <- readland.tps("handaxe.tps", readcurves = TRUE)   #LANDMARK DATA FOR ALL HANDAXE EXAMPLES (SAFE TO IGNORE RESCALING ISSUE)

shape_data_elongated          <- read.csv("elongated.csv", header = TRUE, row.names = 1)  #CLASSIFICATORY DATA FOR ALL ELONGATED EXAMPLES
shape_data_elongated$Artefact <- as.factor(shape_data_elongated$Artefact)
shape_data_tanged             <- read.csv("tanged.csv", header = TRUE, row.names = 1)     #CLASSIFICATORY DATA FOR ALL TANGED EXAMPLES
shape_data_tanged$Artefact    <- as.factor(shape_data_tanged$Artefact)
shape_data_handaxe            <- read.csv("handaxe.csv", header = TRUE, row.names = 1)    #CLASSIFICATORY DATA FOR ALL HANDAXE EXAMPLES
shape_data_handaxe$Artefact   <- as.factor(shape_data_handaxe$Artefact)
shape_data_sliders            <- read.csv("curveslide.csv", header = TRUE)                #SLIDER FILE (PRODUCED USING THE GEOMORPH PACKAGE)

metric_data <- read.csv("measurement_data.csv", header = TRUE, row.names = 1) #METRIC AND CLASSIFICATORY DATA FOR ALL  EXAMPLES
metric_data$Artefact <- as.factor(metric_data$Artefact)

digitisation_error_landmarks  <- readland.tps("digitisation_error.tps", readcurves = TRUE) #LANDMARK DATA FOR A RANDOM EXAMPLE (5 REPLICATES)
digitisation_error_landmarks_data <- read.csv("digitisation_error.csv", header = TRUE, row.names = 1) #CLASSIFICATION: ILLUSTRATION ERROR

digitisation_error_metrics  <- read.csv("measurement_error.csv", header = TRUE, row.names = 1) #METRIC  DATA FOR A RANDOM EXAMPLE (5 REPLICATES)
digitisation_error_metrics$Attempt <- as.factor(digitisation_error_metrics$Attempt)


### STAGE 2: MEASURING INTRA-OBSERVOR ERROR ###

### STAGE 2A: DIGITISATION/LANDMARK ERROR ###

gpa_digi_error <- gpagen(digitisation_error_landmarks, Proj = TRUE, curves =  shape_data_sliders, ProcD = TRUE, surfaces = NULL) #GPA: DIGITISATION ERROR
gpa_digi_error
plot(gpa_digi_error)

gpa_digi_error_df <- geomorph.data.frame(gpa_digi_error, attempt = digitisation_error_landmarks_data$Attempt)
gpaprocD <- procD.lm(coords ~ attempt, data = gpa_digi_error_df) #ANOVA (SHAPE VS INDIVIDUAL)
summary(gpaprocD) #SUMMARY
gpaprocD$aov.table$SS[1]/gpaprocD$aov.table$SS[3]*100 #ERROR AS A PERCENTAGE (8.599812%)

### STAGE 2B: MEASUREMENT ERROR ###

head(digitisation_error_metrics)
statsl  <- describe(digitisation_error_metrics$Length_mm) #DESCRIPTIVE STATISTICS
statsw  <- describe(digitisation_error_metrics$Width_mm) #DESCRIPTIVE STATISTICS
statssf <- describe(digitisation_error_metrics$Scale_Factor) #DESCRIPTIVE STATISTICS

(statsl$se/statsl$mean)   * 100 #FRACTIONAL UNCERTAINITY (LENGTH)
(statsw$se/statsw$mean)   * 100 #FRACTIONAL UNCERTAINITY (WIDTH)
(statssf$se/statssf$mean) * 100 #FRACTIONAL UNCERTAINITY (SCALE FACTOR)


### STAGE 3: GENERALISED PROCRUSTES ANALYSIS (GPA) ###

gpa_elongated <- gpagen(landmarks_elongated, Proj = TRUE, ProcD = TRUE, curves = shape_data_sliders, surfaces = NULL) #GPA: ELONGATED (PROCRUSTES DISTANCE CRITERION)
gpa_elongated$coords #PROCRUSTES COORDINATES FOR ELONGATED EXAMPLES
plot(gpa_elongated) #PLOT: PROCRUSTES COORDINATES ELONGATED

gpa_tanged    <- gpagen(landmarks_tanged, Proj = TRUE, ProcD = TRUE, curves = shape_data_sliders, surfaces = NULL) #GPA: TANGED (PROCRUSTES DISTANCE CRITERION)
gpa_tanged$coords #PROCRUSTES COORDINATES FOR TANGED EXAMPLES
plot(gpa_tanged) #PLOT: PROCRUSTES COORDINATES (TANGED)

gpa_handaxe   <- gpagen(landmarks_handaxe, Proj = TRUE, ProcD = TRUE, curves = shape_data_sliders, surfaces = NULL) #GPA: HANDAXE (PROCRUSTES DISTANCE CRITERION)
gpa_handaxe$coords #PROCRUSTES COORDINATES FOR HANDAXE EXAMPLES
plot(gpa_handaxe) #PLOT: PROCRUSTES COORDINATES (HANDAXES)

### STAGE 4: EXPLORATORY ANALYSIS OF SHAPE DATA (PCA) ###

pca_elongated <- plotTangentSpace(gpa_elongated$coords, cex = 2, warpgrids = FALSE, groups = shape_data_tanged$Artefact, verbose = TRUE) #PCA FOR PROCRUSTES COORDINATES (ELONGATED EXAMPLES: GEOMORPH LAYOUT)
pca_elongated$pc.summary
elongated_ds <- cbind(shape_data_elongated, pca_elongated$pc.scores) #GGPLOT COMPATIBILITY

pca_tanged <- plotTangentSpace(gpa_tanged$coords, warpgrids = FALSE, groups = shape_data_tanged$Artefact, verbose = TRUE) #PCA FOR PROCRUSTES COORDINATES (TANGED: GEOMORPH LAYOUT)
pca_tanged$pc.summary
tanged_ds <- cbind(shape_data_tanged, pca_tanged$pc.scores) #GGPLOT COMPATIBILITY

pca_handaxe <- plotTangentSpace(gpa_handaxe$coords, warpgrids = FALSE, groups = shape_data_handaxe$Artefact, verbose = TRUE) #PCA FOR PROCRUSTES COORDINATES (HANDAXE: GEOMORPH LAYOUT)
pca_handaxe$pc.summary
handaxe_ds <- cbind(shape_data_handaxe, pca_handaxe$pc.scores) #GGPLOT COMPATIBILITY

figure_5 <- ggplot(data = elongated_ds) + geom_point(mapping = aes(x = PC1, y = PC2, colour = Artefact, shape = Class), size = 2) + labs(x = "Principal Component 1 (45.816%)", y = "Principal Component 2 (24.298%)", shape = "Method/skill", colour = "Artefact") + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")) +  scale_shape_manual(values=c(17,15,3,16)) + theme_minimal() 
figure_6 <- ggplot(data = tanged_ds) + geom_point(mapping = aes(x = PC1, y = PC2, colour = Artefact, shape = Class), size = 2) + labs(x = "Principal Component 1 (50.805%)", y = "Principal Component 2 (33.926%)", shape = "Method/skill", colour = "Artefact") + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")) +  scale_shape_manual(values=c(17,15,3,16)) + theme_minimal()
figure_7 <- ggplot(data = handaxe_ds) + geom_point(mapping = aes(x = PC1, y = PC2, colour = Artefact, shape = Class), size = 2) + labs(x = "Principal Component 1 (84.840%)", y = "Principal Component 2 (6.987%)", shape = "Method/skill", colour = "Artefact") + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")) +  scale_shape_manual(values=c(17,15,3,16)) + theme_minimal()

ggsave("Fig5.tiff", plot = figure_5, dpi = 400, units = "mm", height = 150, width = 175)
ggsave("Fig6.tiff", plot = figure_6, dpi = 400, units = "mm", height = 150, width = 175)
ggsave("Fig7.tiff", plot = figure_7, dpi = 400, units = "mm", height = 150, width = 175)


### STAGE 5: ANALYTICAL FRAMEWORK FOR THE SHAPE DATA ###

### STAGE 5A: ELONGATED EXAMPLES ###

elongated1 <- gpa_elongated$coords[, , 1:18] #PROCRUSTES COORDINATES: FIRST 18 (ARTEFACT 1)
df_elongated1 <- geomorph.data.frame(shape = elongated1, class = shape_data_elongated$Class[1:18], artefact = shape_data_elongated$Artefact[1:18]) #DATA FRAME CREATION
df_elongated1 #SUMMARY: DATA FRAME (GEOMORPH FORMAT)
E1 <- procD.lm(shape ~ class, data = df_elongated1) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 1)
summary(E1) #ANOVA SUMMARY

elongated1pcs <- as.data.frame(pca_elongated$pc.scores[1:18, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
elongated1class <- shape_data_elongated$Class[1:18]
elongated1pcs <- cbind(elongated1pcs, elongated1class)
elongated1pcs <- rename(elongated1pcs, Class = elongated1class)
elongated1lda <- lda(Class ~ ., data = elongated1pcs)
elongated1ldapredict <- predict(elongated1lda)
elongated1ldaplot <- cbind(elongated1ldapredict$x[, 1:3], elongated1pcs)
elongated1ldaplotggplot <- ggplot(elongated1ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.title.y = element_text(size = 6), legend.position = "none",  axis.title.x = element_text(size = 6))

elongated2 <- gpa_elongated$coords[, , 19:36] #PROCRUSTES COORDINATES: SECOND 18 (ARTEFACT 2)
df_elongated2 <- geomorph.data.frame(shape = elongated2, class = shape_data_elongated$Class[19:36], artefact = shape_data_elongated$Artefact[19:36]) #DATA FRAME CREATION
df_elongated2 #SUMMARY: DATA FRAME
E2 <- procD.lm(shape ~ class, data = df_elongated2) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 2)
summary(E2) #ANOVA SUMMARY
 
elongated2pcs <- as.data.frame(pca_elongated$pc.scores[19:36, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
elongated2class <- shape_data_elongated$Class[19:36]
elongated2pcs <- cbind(elongated2pcs, elongated2class)
elongated2pcs <- rename(elongated2pcs, Class = elongated2class)
elongated2lda <- lda(Class ~ ., elongated2pcs)
elongated2ldapredict <- predict(elongated2lda)
elongated2ldaplot <- cbind(elongated2ldapredict$x[, 1:3], elongated2pcs)
elongated2ldaplotggplot <- ggplot(elongated2ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

elongated3 <- gpa_elongated$coords[, , 37:54] #PROCRUSTES COORDINATES: THIRD 18 (ARTEFACT 3)
df_elongated3 <- geomorph.data.frame(shape = elongated3, class = shape_data_elongated$Class[37:54], artefact = shape_data_elongated$Artefact[37:54]) #DATA FRAME CREATION
df_elongated3 #SUMMARY: DATA FRAME
E3 <- procD.lm(shape ~ class, data = df_elongated3) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 3)
summary(E3) #ANOVA SUMMARY

elongated3pcs <- as.data.frame(pca_elongated$pc.scores[37:54, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
elongated3class <- shape_data_elongated$Class[37:54]
elongated3pcs <- cbind(elongated3pcs, elongated3class)
elongated3pcs <- rename(elongated3pcs, Class = elongated3class)
elongated3lda <- lda(Class ~ ., elongated3pcs)
elongated3ldapredict <- predict(elongated3lda)
elongated3ldaplot <- cbind(elongated3ldapredict$x[, 1:3], elongated3pcs)
elongated3ldaplotggplot <- ggplot(elongated3ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

elongated4 <- gpa_elongated$coords[, , 55:72] #PROCRUSTES COORDINATES: FOURTH 18 (ARTEFACT 4)
df_elongated4 <- geomorph.data.frame(shape = elongated4, class = shape_data_elongated$Class[55:72], artefact = shape_data_elongated$Artefact[55:72]) #DATA FRAME CREATION
df_elongated4 #SUMMARY: DATA FRAME
E4 <- procD.lm(shape ~ class, data = df_elongated4) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 4)
summary(E4) #ANOVA SUMMARY

elongated4pcs <- as.data.frame(pca_elongated$pc.scores[55:72, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
elongated4class <- shape_data_elongated$Class[55:72]
elongated4pcs <- cbind(elongated4pcs, elongated4class)
elongated4pcs <- rename(elongated4pcs, Class = elongated4class)
elongated4lda <- lda(Class ~ ., elongated4pcs)
elongated4ldapredict <- predict(elongated4lda)
elongated4ldaplot <- cbind(elongated4ldapredict$x[, 1:3], elongated4pcs)
elongated4ldaplotggplot <- ggplot(elongated4ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

elongated5 <- gpa_elongated$coords[, , 73:90] #PROCRUSTES COORDINATES: FIFTH 18 (ARTEFACT 5)
df_elongated5 <- geomorph.data.frame(shape = elongated5, class = shape_data_elongated$Class[73:90], artefact = shape_data_elongated$Artefact[73:90]) #DATA FRAME CREATION
df_elongated5 #SUMMARY: DATA FRAME
E5 <- procD.lm(shape ~ class, data = df_elongated5) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 5)
summary(E5) #ANOVA SUMMARY

elongated5pcs <- as.data.frame(pca_elongated$pc.scores[73:90, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
elongated5class <- shape_data_elongated$Class[73:90]
elongated5pcs <- cbind(elongated5pcs, elongated5class)
elongated5pcs <- rename(elongated5pcs, Class = elongated5class)
elongated5lda <- lda(Class ~ ., elongated5pcs)
elongated5ldapredict <- predict(elongated5lda)
elongated5ldaplot <- cbind(elongated5ldapredict$x[, 1:3], elongated5pcs)
elongated5ldaplotggplot <- ggplot(elongated5ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

### STAGE 5B: TANGED EXAMPLES ###

tanged1 <- gpa_tanged$coords[, , 1:18] #PROCRUSTES COORDINATES: FIRST 18 (ARTEFACT 1)
df_tanged1 <- geomorph.data.frame(shape = tanged1, class = shape_data_tanged$Class[1:18], artefact = shape_data_tanged$Artefact[1:18]) #DATA FRAME CREATION
df_tanged1 #SUMMARY: DATA FRAME
T1 <- procD.lm(shape ~ class, data = df_tanged1) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 1)
summary(T1) #ANOVA SUMMARY

tanged1pcs <- as.data.frame(pca_tanged$pc.scores[1:18, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
tanged1class <- shape_data_tanged$Class[1:18]
tanged1pcs <- cbind(tanged1pcs, tanged1class)
tanged1pcs <- rename(tanged1pcs, Class = tanged1class)
tanged1lda <- lda(Class ~ ., tanged1pcs)
tanged1ldapredict <- predict(tanged1lda)
tanged1ldaplot <- cbind(tanged1ldapredict$x[, 1:3], tanged1pcs)
tanged1ldaplotggplot <- ggplot(tanged1ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

tanged2 <- gpa_tanged$coords[, , 19:36] #PROCRUSTES COORDINATES: SECOND 18 (ARTEFACT 2)
df_tanged2 <- geomorph.data.frame(shape = tanged2, class = shape_data_tanged$Class[19:36], artefact = shape_data_tanged$Artefact[19:36]) #DATA FRAME CREATION
df_tanged2 #SUMMARY: DATA FRAME
T2 <- procD.lm(shape ~ class, data = df_tanged2) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 2)
summary(T2) #ANOVA SUMMARY

tanged2pcs <- as.data.frame(pca_tanged$pc.scores[19:36, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
tanged2class <- shape_data_tanged$Class[19:36]
tanged2pcs <- cbind(tanged2pcs, tanged2class)
tanged2pcs <- rename(tanged2pcs, Class = tanged2class)
tanged2lda <- lda(Class ~ ., tanged2pcs)
tanged2ldapredict <- predict(tanged2lda)
tanged2ldaplot <- cbind(tanged2ldapredict$x[, 1:3], tanged2pcs)
tanged2ldaplotggplot <- ggplot(tanged2ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

tanged3 <- gpa_tanged$coords[, , 37:54] #PROCRUSTES COORDINATES: THIRD 18 (ARTEFACT 3)
df_tanged3 <- geomorph.data.frame(shape = tanged3, class = shape_data_tanged$Class[37:54], artefact = shape_data_tanged$Artefact[37:54]) #DATA FRAME CREATION
df_tanged3 #SUMMARY: DATA FRAME
T3 <- procD.lm(shape ~ class, data = df_tanged3) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 3)
summary(T3) #ANOVA SUMMARY

tanged3pcs <- as.data.frame(pca_tanged$pc.scores[37:54, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
tanged3class <- shape_data_tanged$Class[37:54]
tanged3pcs <- cbind(tanged3pcs, tanged3class)
tanged3pcs <- rename(tanged3pcs, Class = tanged3class)
tanged3lda <- lda(Class ~ ., tanged3pcs)
tanged3ldapredict <- predict(tanged3lda)
tanged3ldaplot <- cbind(tanged3ldapredict$x[, 1:3], tanged3pcs)
tanged3ldaplotggplot <- ggplot(tanged3ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

tanged4 <- gpa_tanged$coords[, , 55:72] #PROCRUSTES COORDINATES: FOURTH 18 (ARTEFACT 4)
df_tanged4 <- geomorph.data.frame(shape = tanged4, class = shape_data_tanged$Class[55:72], artefact = shape_data_tanged$Artefact[55:72]) #DATA FRAME CREATION
df_tanged4 #SUMMARY: DATA FRAME
T4 <- procD.lm(shape ~ class, data = df_tanged4) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 4)
summary(T4) #ANOVA SUMMARY

tanged4pcs <- as.data.frame(pca_tanged$pc.scores[55:72, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
tanged4class <- shape_data_tanged$Class[55:72]
tanged4pcs <- cbind(tanged4pcs, tanged4class)
tanged4pcs <- rename(tanged4pcs, Class = tanged4class)
tanged4lda <- lda(Class ~ ., tanged4pcs)
tanged4ldapredict <- predict(tanged4lda)
tanged4ldaplot <- cbind(tanged4ldapredict$x[, 1:3], tanged4pcs)
tanged4ldaplotggplot <- ggplot(tanged4ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

tanged5 <- gpa_tanged$coords[, , 73:90] #PROCRUSTES COORDINATES: FIFTH 18 (ARTEFACT 5)
df_tanged5 <- geomorph.data.frame(shape = tanged5, class = shape_data_tanged$Class[73:90], artefact = shape_data_tanged$Artefact[73:90]) #DATA FRAME CREATION
df_tanged5 #SUMMARY: DATA FRAME
T5 <- procD.lm(shape ~ class, data = df_tanged5) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 5)
summary(T5) #ANOVA SUMMARY

tanged5pcs <- as.data.frame(pca_tanged$pc.scores[73:90, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
tanged5class <- shape_data_tanged$Class[73:90]
tanged5pcs <- cbind(tanged5pcs, tanged5class)
tanged5pcs <- rename(tanged5pcs, Class = tanged5class)
tanged5lda <- lda(Class ~ ., tanged5pcs)
tanged5ldapredict <- predict(tanged5lda)
tanged5ldaplot <- cbind(tanged5ldapredict$x[, 1:3], tanged5pcs)
tanged5ldaplotggplot <- ggplot(tanged5ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

### STAGE 5C: HANDAXE EXAMPLES ###

handaxe1 <- gpa_handaxe$coords[, , 1:18] #PROCRUSTES COORDINATES: FIRST 18 (ARTEFACT 1)
df_handaxe1 <- geomorph.data.frame(shape = handaxe1, class = shape_data_handaxe$Class[1:18], artefact = shape_data_handaxe$Artefact[1:18]) #DATA FRAME CREATION
df_handaxe1 #SUMMARY: DATA FRAME
H1 <- procD.lm(shape ~ class, data = df_handaxe1) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 1)
summary(H1) #ANOVA SUMMARY

handaxe1pcs <- as.data.frame(pca_handaxe$pc.scores[1:18, 1:10]) #91% CUMULATIVE SHAPE VARIANCE
handaxe1class <- shape_data_handaxe$Class[1:18]
handaxe1pcs <- cbind(handaxe1pcs, handaxe1class)
handaxe1pcs <- rename(handaxe1pcs, Class = handaxe1class)
handaxe1lda <- lda(Class ~ ., handaxe1pcs)
handaxe1ldapredict <- predict(handaxe1lda)
handaxe1ldaplot <- cbind(handaxe1ldapredict$x[, 1:3], handaxe1pcs)
handaxe1ldaplotggplot <- ggplot(handaxe1ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

handaxe2 <- gpa_handaxe$coords[, , 19:36] #PROCRUSTES COORDINATES: SECOND 18 (ARTEFACT 2)
df_handaxe2 <- geomorph.data.frame(shape = handaxe2, class = shape_data_handaxe$Class[19:36], artefact = shape_data_handaxe$Artefact[19:36]) #DATA FRAME CREATION
df_handaxe2 #SUMMARY: DATA FRAME
H2 <- procD.lm(shape ~ class, data = df_handaxe2) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 2)
summary(H2) #ANOVA SUMMARY

handaxe2pcs <- as.data.frame(pca_handaxe$pc.scores[19:36, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
handaxe2class <- shape_data_handaxe$Class[19:36]
handaxe2pcs <- cbind(handaxe2pcs, handaxe2class)
handaxe2pcs <- rename(handaxe2pcs, Class = handaxe2class)
handaxe2lda <- lda(Class ~ ., handaxe2pcs)
handaxe2ldapredict <- predict(handaxe2lda)
handaxe2ldaplot <- cbind(handaxe2ldapredict$x[, 1:3], handaxe2pcs)
handaxe2ldaplotggplot <- ggplot(handaxe2ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

handaxe3 <- gpa_handaxe$coords[, , 37:54] #PROCRUSTES COORDINATES: THIRD 18 (ARTEFACT 3)
df_handaxe3 <- geomorph.data.frame(shape = handaxe3, class = shape_data_handaxe$Class[37:54], artefact = shape_data_handaxe$Artefact[37:54]) #DATA FRAME CREATION
df_handaxe3 #SUMMARY: DATA FRAME
H3 <- procD.lm(shape ~ class, data = df_handaxe3) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 3)
summary(H3) #ANOVA SUMMARY

handaxe3pcs <- as.data.frame(pca_handaxe$pc.scores[37:54, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
handaxe3class <- shape_data_handaxe$Class[37:54]
handaxe3pcs <- cbind(handaxe3pcs, handaxe3class)
handaxe3pcs <- rename(handaxe3pcs, Class = handaxe3class)
handaxe3lda <- lda(Class ~ ., handaxe3pcs)
handaxe3ldapredict <- predict(handaxe3lda)
handaxe3ldaplot <- cbind(handaxe3ldapredict$x[, 1:3], handaxe3pcs)
handaxe3ldaplotggplot <- ggplot(handaxe3ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

handaxe4 <- gpa_handaxe$coords[, , 55:72] #PROCRUSTES COORDINATES: FOURTH 18 (ARTEFACT 4)
df_handaxe4 <- geomorph.data.frame(shape = handaxe4, class = shape_data_handaxe$Class[55:72], artefact = shape_data_handaxe$Artefact[55:72]) #DATA FRAME CREATION
df_handaxe4 #SUMMARY: DATA FRAME
H4 <- procD.lm(shape ~ class, data = df_handaxe4) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 4)
summary(H4) #ANOVA SUMMARY

handaxe4pcs <- as.data.frame(pca_handaxe$pc.scores[55:72, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
handaxe4class <- shape_data_handaxe$Class[55:72]
handaxe4pcs <- cbind(handaxe4pcs, handaxe4class)
handaxe4pcs <- rename(handaxe4pcs, Class = handaxe4class)
handaxe4lda <- lda(Class ~ ., handaxe4pcs)
handaxe4ldapredict <- predict(handaxe4lda)
handaxe4ldaplot <- cbind(handaxe4ldapredict$x[, 1:3], handaxe4pcs)
handaxe4ldaplotggplot <- ggplot(handaxe4ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

handaxe5 <- gpa_handaxe$coords[, , 73:90] #PROCRUSTES COORDINATES: FIFTH 18 (ARTEFACT 5)
df_handaxe5 <- geomorph.data.frame(shape = handaxe5, class = shape_data_handaxe$Class[73:90], artefact = shape_data_handaxe$Artefact[73:90]) #DATA FRAME CREATION
df_handaxe5 #SUMMARY: DATA FRAME
H5 <- procD.lm(shape ~ class, data = df_handaxe5) #PROCRUSTES ANOVA: SHAPE VS CLASS (ARTEFACT 5)
summary(H5) #ANOVA SUMMARY

handaxe5pcs <- as.data.frame(pca_handaxe$pc.scores[73:90, 1:10]) #95% CUMULATIVE SHAPE VARIANCE
handaxe5class <- shape_data_handaxe$Class[73:90]
handaxe5pcs <- cbind(handaxe5pcs, handaxe5class)
handaxe5pcs <- rename(handaxe5pcs, Class = handaxe5class)
handaxe5lda <- lda(Class ~ ., handaxe5pcs)
handaxe5ldapredict <- predict(handaxe5lda)
handaxe5ldaplot <- cbind(handaxe5ldapredict$x[, 1:3], handaxe5pcs)
handaxe5ldaplotggplot <- ggplot(handaxe5ldaplot, aes(LD1, LD2)) + geom_point(aes(shape = Class), size = 2) + labs(x = "LDA 1", y = "LDA 2") + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",  axis.title.y = element_text(size = 6),  axis.title.x = element_text(size = 6))

figure_8 <- plot_grid(elongated1ldaplotggplot, elongated2ldaplotggplot, elongated3ldaplotggplot, elongated4ldaplotggplot, elongated5ldaplotggplot, tanged1ldaplotggplot, tanged2ldaplotggplot, tanged3ldaplotggplot, tanged4ldaplotggplot, tanged5ldaplotggplot, handaxe1ldaplotggplot, handaxe2ldaplotggplot, handaxe3ldaplotggplot, handaxe4ldaplotggplot, handaxe5ldaplotggplot, ncol = 5, align = 'v') #synthesis of the four figures
plot(figure_8)
annotate_figure(figure_8, left = text_grob("     Handaxe Illustrations                 Tanged Illustrations                 Elongated Illustrations", rot = 90, size = 10))
ggsave("Fig8.tiff", plot = ggplot2::last_plot(), units = "mm", dpi = 300, height = 150, width = 250)

### STAGE 6: EXPLORATORY FRAMEWORK FOR THE MEASUREMENT DATA ###

metric_data_elongated <- metric_data[which(metric_data$Type=="Elongated"),]
metric_data_handaxe   <- metric_data[which(metric_data$Type=="Handaxe"),]
metric_data_tanged    <- metric_data[which(metric_data$Type=="Tanged"),]

figure_9a <- ggplot(metric_data_elongated, aes(Length_mm, Width_mm, colour = Artefact, shape = Class)) + geom_point() + facet_grid(cols = vars(Artefact), scales = "free", labeller=label_both) + labs(x = "Length (mm)", y = "Width (mm)") + theme_minimal() + scale_shape_manual(values=c(17,15,3,16)) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
figure_9b <- ggplot(metric_data_tanged, aes(Length_mm, Width_mm, colour = Artefact, shape = Class)) + geom_point() + facet_grid(cols = vars(Artefact), scales = "free", labeller=label_both) + labs(x = "Length (mm)", y = "Width (mm)") +  scale_shape_manual(values=c(17,15,3,16)) + theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
figure_9c <- ggplot(metric_data_handaxe, aes(Length_mm, Width_mm, colour = Artefact, shape = Class)) + geom_point() + facet_grid(cols = vars(Artefact), scales = "free", labeller=label_both) + labs(x = "Length (mm)", y = "Width (mm)") +  scale_shape_manual(values=c(17,15,3,16)) +  theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
figure_9 <- plot_grid(figure_9a, figure_9b, figure_9c, labels= "AUTO", ncol = 1, align = 'v') #synthesis of the four figures
plot(figure_9)
ggsave("Fig9.tiff", plot = last_plot(), dpi = 400, units = "mm", height = 250, width = 200)

### STAGE 7: ANALYTICAL FRAMEWORK FOR THE MEASUREMENT DATA ###

### STAGE 7A: ELONGATED EXAMPLES ###

metric_data_elongated <- arrange(metric_data_elongated, Artefact, Class)

metric_data_elongated[1:18,]  %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[1:18,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[19:36,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[19:36,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[37:54,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[37:54,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[55:72,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[55:72,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[73:90,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[73:90,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)

metric_data_elongated[1:18,]  %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[1:18,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[19:36,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[19:36,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[37:54,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[37:54,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[55:72,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[55:72,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[73:90,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_elongated[73:90,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)

shapiro.test(metric_data_elongated$Length_mm[2:4]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (ELONGATED: ARTEFACT 1)
shapiro.test(metric_data_elongated$Width_mm[2:4]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (ELONGATED: ARTEFACT 1)
shapiro.test(metric_data_elongated$Length_mm[5:17]) #NORMALITY TEST FOR LENGTH: NOVICE (ELONGATED: ARTEFACT 1)
shapiro.test(metric_data_elongated$Width_mm[5:17]) #NORMALITY TEST FOR WIDTH: NOVICE (ELONGATED: ARTEFACT 1)
shapiro.test(metric_data_elongated$Length_mm[20:22]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (ELONGATED: ARTEFACT 2)
shapiro.test(metric_data_elongated$Width_mm[20:22]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (ELONGATED: ARTEFACT 2)
shapiro.test(metric_data_elongated$Length_mm[23:35]) #NORMALITY TEST FOR LENGTH: NOVICE (ELONGATED: ARTEFACT 2)
shapiro.test(metric_data_elongated$Width_mm[23:35]) #NORMALITY TEST FOR WIDTH: NOVICE (ELONGATED: ARTEFACT 2)
shapiro.test(metric_data_elongated$Length_mm[38:40]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (ELONGATED: ARTEFACT 3)
shapiro.test(metric_data_elongated$Width_mm[38:40]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (ELONGATED: ARTEFACT 3)
shapiro.test(metric_data_elongated$Length_mm[41:53]) #NORMALITY TEST FOR LENGTH: NOVICE (ELONGATED: ARTEFACT 3)
shapiro.test(metric_data_elongated$Width_mm[41:53]) #NORMALITY TEST FOR WIDTH: NOVICE (ELONGATED: ARTEFACT 3)
shapiro.test(metric_data_elongated$Length_mm[56:58]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (ELONGATED: ARTEFACT 4)
shapiro.test(metric_data_elongated$Width_mm[56:58]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (ELONGATED: ARTEFACT 4)
shapiro.test(metric_data_elongated$Length_mm[59:71]) #NORMALITY TEST FOR LENGTH: NOVICE (ELONGATED: ARTEFACT 4)
shapiro.test(metric_data_elongated$Width_mm[59:71]) #NORMALITY TEST FOR WIDTH: NOVICE (ELONGATED: ARTEFACT 4)
shapiro.test(metric_data_elongated$Length_mm[74:76]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (ELONGATED: ARTEFACT 5)
shapiro.test(metric_data_elongated$Width_mm[74:76]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (ELONGATED: ARTEFACT 5)
shapiro.test(metric_data_elongated$Length_mm[77:89]) #NORMALITY TEST FOR LENGTH: NOVICE (ELONGATED: ARTEFACT 5)
shapiro.test(metric_data_elongated$Width_mm[77:89]) #NORMALITY TEST FOR WIDTH: NOVICE (ELONGATED: ARTEFACT 5)

summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[1:18,])) #MANOVA ELONGATED (1ST ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[1:18,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[19:36,]))  #MANOVA ELONGATED (2ND ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[19:36,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[37:54,]))  #MANOVA ELONGATED (3RD ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[37:54,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[55:72,]))  #MANOVA ELONGATED (4TH ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[55:72,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[73:90,]))  #MANOVA ELONGATED (5TH ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_elongated[73:90,]))  #INDIVIDUAL RESPONSES

### STAGE 7B: HANDAXE EXAMPLES ###

metric_data_handaxe <- arrange(metric_data_handaxe, Artefact, Class)

metric_data_handaxe[1:18,]  %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[1:18,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[19:36,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[19:36,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[37:54,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[37:54,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[55:72,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[55:72,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[73:90,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[73:90,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)

metric_data_handaxe[1:18,]  %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[1:18,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[19:36,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[19:36,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[37:54,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[37:54,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[55:72,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[55:72,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[73:90,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_handaxe[73:90,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)

shapiro.test(metric_data_handaxe$Length_mm[2:4]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (HANDAXE: ARTEFACT 1)
shapiro.test(metric_data_handaxe$Width_mm[2:4]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (HANDAXE: ARTEFACT 1)
shapiro.test(metric_data_handaxe$Length_mm[5:17]) #NORMALITY TEST FOR LENGTH: NOVICE (HANDAXE: ARTEFACT 1)
shapiro.test(metric_data_handaxe$Width_mm[5:17]) #NORMALITY TEST FOR WIDTH: NOVICE (HANDAXE: ARTEFACT 1)
shapiro.test(metric_data_handaxe$Length_mm[20:22]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (HANDAXE: ARTEFACT 2)
shapiro.test(metric_data_handaxe$Width_mm[20:22]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (HANDAXE: ARTEFACT 2)
shapiro.test(metric_data_handaxe$Length_mm[23:35]) #NORMALITY TEST FOR LENGTH: NOVICE (HANDAXE: ARTEFACT 2)
shapiro.test(metric_data_handaxe$Width_mm[23:35]) #NORMALITY TEST FOR WIDTH: NOVICE (HANDAXE: ARTEFACT 2)
shapiro.test(metric_data_handaxe$Length_mm[38:40]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (HANDAXE: ARTEFACT 3)
shapiro.test(metric_data_handaxe$Width_mm[38:40]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (HANDAXE: ARTEFACT 3)
shapiro.test(metric_data_handaxe$Length_mm[41:53]) #NORMALITY TEST FOR LENGTH: NOVICE (HANDAXE: ARTEFACT 3)
shapiro.test(metric_data_handaxe$Width_mm[41:53]) #NORMALITY TEST FOR WIDTH: NOVICE (HANDAXE: ARTEFACT 3) - NON-NORMAL RESULT
shapiro.test(metric_data_handaxe$Length_mm[56:58]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (HANDAXE: ARTEFACT 4)
shapiro.test(metric_data_handaxe$Width_mm[56:58]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (HANDAXE: ARTEFACT 4)
shapiro.test(metric_data_handaxe$Length_mm[59:71]) #NORMALITY TEST FOR LENGTH: NOVICE (HANDAXE: ARTEFACT 4)
shapiro.test(metric_data_handaxe$Width_mm[59:71]) #NORMALITY TEST FOR WIDTH: NOVICE (HANDAXE: ARTEFACT 4)
shapiro.test(metric_data_handaxe$Length_mm[74:76]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (HANDAXE: ARTEFACT 5)
shapiro.test(metric_data_handaxe$Width_mm[74:76]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (HANDAXE: ARTEFACT 5)
shapiro.test(metric_data_handaxe$Length_mm[77:89]) #NORMALITY TEST FOR LENGTH: NOVICE (HANDAXE: ARTEFACT 5)
shapiro.test(metric_data_handaxe$Width_mm[77:89]) #NORMALITY TEST FOR WIDTH: NOVICE (handaxe: ARTEFACT 5)

summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[1:18,])) #MANOVA HANDAXE (1ST ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[1:18,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[19:36,]))  #MANOVA HANDAXE (2ND ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[19:36,]))  #INDIVIDUAL RESPONSES
adonis(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[37:54,])  #NON-PARAMETRIC MANOVA HANDAXE (3RD ARTEFACT)
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[55:72,]))   #MANOVA HANDAXE (4TH ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[55:72,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[73:90,]))   #MANOVA HANDAXE (5TH ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_handaxe[73:90,]))  #INDIVIDUAL RESPONSES


### STAGE 7C: TANGED EXAMPLES ###

metric_data_tanged <- arrange(metric_data_tanged, Artefact, Class)

metric_data_tanged[1:18,]  %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[1:18,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[19:36,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[19:36,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[37:54,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[37:54,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[55:72,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[55:72,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[73:90,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Length_mm, na.rm = TRUE), sd = sd(Length_mm, na.rm = TRUE), min = min(Length_mm, na.rm = TRUE), max = max(Length_mm, na.rm = TRUE), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[73:90,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Length_mm, na.rm = TRUE)/mean(Length_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)

metric_data_tanged[1:18,]  %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[1:18,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[19:36,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[19:36,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[37:54,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[37:54,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[55:72,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[55:72,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[73:90,] %>% group_by(Class) %>% summarise(count = n(), mean = mean(Width_mm, na.rm = TRUE), sd = sd(Width_mm, na.rm = TRUE), min = min(Width_mm, na.rm = TRUE), max = max(Width_mm, na.rm = TRUE), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)
metric_data_tanged[73:90,]  %>% filter(Class=="Novice" | Class == "Intermediate" | Class == "Professional") %>% summarise(count = n(), cv = sd(Width_mm, na.rm = TRUE)/mean(Width_mm, na.rm = TRUE)*100) %>% dplyr::mutate_if(is.numeric, format, 1)

shapiro.test(metric_data_tanged$Length_mm[2:4]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (TANGED: ARTEFACT 1)
shapiro.test(metric_data_tanged$Width_mm[2:4]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (TANGED: ARTEFACT 1)
shapiro.test(metric_data_tanged$Length_mm[5:17]) #NORMALITY TEST FOR LENGTH: NOVICE (TANGED: ARTEFACT 1)
shapiro.test(metric_data_tanged$Width_mm[5:17]) #NORMALITY TEST FOR WIDTH: NOVICE (TANGED: ARTEFACT 1)
shapiro.test(metric_data_tanged$Length_mm[20:22]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (TANGED: ARTEFACT 2)
shapiro.test(metric_data_tanged$Width_mm[20:22]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (TANGED: ARTEFACT 2)
shapiro.test(metric_data_tanged$Length_mm[23:35]) #NORMALITY TEST FOR LENGTH: NOVICE (TANGED: ARTEFACT 2)
shapiro.test(metric_data_tanged$Width_mm[23:35]) #NORMALITY TEST FOR WIDTH: NOVICE (TANGED: ARTEFACT 2)
shapiro.test(metric_data_tanged$Length_mm[38:40]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (TANGED: ARTEFACT 3)
shapiro.test(metric_data_tanged$Width_mm[38:40]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (TANGED: ARTEFACT 3)
shapiro.test(metric_data_tanged$Length_mm[41:53]) #NORMALITY TEST FOR LENGTH: NOVICE (TANGED: ARTEFACT 3)
shapiro.test(metric_data_tanged$Width_mm[41:53]) #NORMALITY TEST FOR WIDTH: NOVICE (TANGED: ARTEFACT 3)
shapiro.test(metric_data_tanged$Length_mm[56:58]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (TANGED: ARTEFACT 4)
shapiro.test(metric_data_tanged$Width_mm[56:58]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (TANGED: ARTEFACT 4)
shapiro.test(metric_data_tanged$Length_mm[59:71]) #NORMALITY TEST FOR LENGTH: NOVICE (TANGED: ARTEFACT 4)
shapiro.test(metric_data_tanged$Width_mm[59:71]) #NORMALITY TEST FOR WIDTH: NOVICE (TANGED: ARTEFACT 4)
shapiro.test(metric_data_tanged$Length_mm[74:76]) #NORMALITY TEST FOR LENGTH: INTERMEDIATE (TANGED: ARTEFACT 5)
shapiro.test(metric_data_tanged$Width_mm[74:76]) #NORMALITY TEST FOR WIDTH: INTERMEDIATE (TANGED: ARTEFACT 5)
shapiro.test(metric_data_tanged$Length_mm[77:89]) #NORMALITY TEST FOR LENGTH: NOVICE (TANGED: ARTEFACT 5)
shapiro.test(metric_data_tanged$Width_mm[77:89]) #NORMALITY TEST FOR WIDTH: NOVICE (TANGED: ARTEFACT 5)

summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[1:18,])) #MANOVA TANGED (1ST ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[1:18,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[19:36,]))  #MANOVA TANGED (2ND ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[19:36,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[37:54,]))  #MANOVA TANGED (3RD ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[37:54,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[55:72,]))  #MANOVA TANGED (4TH ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[55:72,]))  #INDIVIDUAL RESPONSES
summary(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[73:90,]))  #MANOVA TANGED (5TH ARTEFACT)
summary.aov(manova(cbind(Length_mm, Width_mm) ~ Class, data = metric_data_tanged[73:90,]))  #INDIVIDUAL RESPONSES

########
