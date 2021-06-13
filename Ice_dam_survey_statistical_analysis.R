###################################################################################
# Basic GLM to clarify effect of river parameters on ice dam formation probability#
###################################################################################

#Libraries
library(dplyr) #For data wrangling
library(ggplot2) #For plotting data
library(cowplot) #For placing plots in a grid
library(GGally) #For making correlation heatmaps
library(ROCR) #For making the ROC curve
library(MuMIn) #For model subsetting and model selection
library(boot) #For cross-validation

# Setwd
setwd("A:/PhD/Analysis/Ice_dam_survey")

# read data
data = read.csv("WidthSlopeSinu.csv", sep=";")#Importing CSV file, make sure the separator is correct

#Cursory inspection of data
glimpse(data)

#Remove irrelevant columns from data
data <- select(data, -c(name, 
                        height, 
                        north_Y_, 
                        east_X_,
                        layer,
                        Dam_width,
                        Dam_Slope,
                        Dam_SlopeMovWin, 
                        Dam_SinuosityMovWin, 
                        Dam_Sinuosity140,
                        Dam_SinuosityCustom))
#Rename poorly named column
data <- rename(data, "Distance"="ï..dist")

#Omit rows with NA entries
data <-na.omit(data)

#Extract columns with continuous, as opposed to binary data
continuous <- select_if(data, is.numeric)

#Summary statistics of continuous data
summary(continuous)

#Standardize continuous variables
#data <- data %>%
#  mutate_if(is.numeric, funs(as.numeric(scale(.))))
#head(data) #Inspect outcome

#Convert chr columns to factors
data <- data %>% mutate_if(is.character,as.factor)
levels(data$Type) #Inspect factor levels in Type

#Remove rows that are neither "No Dam" nor "Ice Dam"
data <- filter(data, Type!="") %>% droplevels
levels(data$Type) #Inspect factor levels in Type
glimpse(data)

#Plot number of dams and no dams
ggplot(data, aes(Type))+
  geom_bar()

#Plot Ice dam/No Dam vs all continuous variables
graph_box_plot <- lapply(names(continuous),
                function(y)
                  ggplot(data, aes(x=Type, get(y)))+
                  geom_boxplot()+
                  stat_summary(fun=mean,
                               geom="point",
                               size=3,
                               color="steelblue")+
                  ylab(y)+
                  theme_classic())
graph_box_plot

#Plot distributions by dam status
graph_density_by_Type <-lapply(names(continuous),
                               function(x)
                                 ggplot(data, aes(get(x)))+
                                 geom_density(aes(color=Type), alpha=0.5)+
                                 xlab(x)+
                                 theme_classic())
graph_density_by_Type

plot_grid(graph_density_by_Type[[2]],
          graph_density_by_Type[[3]],
          graph_density_by_Type[[4]],
          graph_density_by_Type[[5]],
          graph_density_by_Type[[6]],
          graph_density_by_Type[[7]],
          ncol=2)


#Carry out one-way anova for each of the continuous variables
#This is to see if there is statistically significant difference between variables depending on presence of icedam
anova <- lapply(names(continuous),
                function(x)
                  aov(get(x)~Type, data)
                )
for(i in 1:length(anova)){
  print(names(continuous)[[i]])
  print(summary(anova[[i]]))
}


#Check whether variables depend on each other
graph_all_vs_slope_cleaned <- lapply(names(continuous),
                function(x)
                  ggplot(data, aes(get(x), y=Slope_cleaned))+
                  geom_point(aes(color=Type),
                             size=0.5)+
                  stat_smooth(method='lm',
                              formula=y~poly(x,2),
                              se = TRUE,
                              aes(color=Type))+
                  xlab(x)+
                  theme_classic()
)
graph_all_vs_slope_cleaned

graph_all_vs_width <- lapply(names(continuous),
                                     function(x)
                                       ggplot(data, aes(get(x), y=Width))+
                                       geom_point(aes(color=Type),
                                                  size=0.5)+
                                       stat_smooth(method='lm',
                                                   formula=y~poly(x,2),
                                                   se = TRUE,
                                                   aes(color=Type))+
                                       xlab(x)+
                                       theme_classic()
)
graph_all_vs_width

graph_all_vs_distance <- lapply(names(continuous),
                             function(x)
                               ggplot(data, aes(get(x), y=distance))+
                               geom_point(aes(color=Type),
                                          size=0.5)+
                               stat_smooth(method='lm',
                                           formula=y~poly(x,2),
                                           se = TRUE,
                                           aes(color=Type))+
                               xlab(x)+
                               theme_classic()
)
graph_all_vs_distance

#Check correlation between variables
data_corr <- data %>%
  mutate_if(is.numeric, funs(as.numeric(scale(.))))
corr <- data.frame(lapply(data_corr, as.integer))
ggcorr(corr,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")

#Fit generalized linear models
formula <- Type~.*.
glm_all <- glm(formula, data=data, family='binomial', na.action=na.fail)
summary(glm_all)

data_1 <- data %>%
  select(Type, Width, Slope_cleaned, Sinuosity_MovWin)
glm_1 <- glm(formula, data=data_1, family='binomial', na.action=na.fail)
summary(glm_1)

data_2 <- data %>%
  select(Type, Width, Slope_cleaned, Sinuosity140)
glm_2 <- glm(formula, data=data_2, family='binomial', na.action=na.fail)
summary(glm_2)

data_3 <- data %>%
  select(Type, Width, Slope_cleaned)
glm_3 <- glm(formula, data=data_3, family='binomial', na.action=na.fail)
summary(glm_3)

dredgeing <- dredge(glm_1, beta="none")
dredgeing
model.avg(dredgeing, beta="none")
best_model <- glm_3
data_best <- data_3

#Assess performance of model using a confusion matrix
predict <- predict(glm_1,data_best, type='response')
table_mat <- table(data_best$Type, predict>0.5)
table_mat
accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test
precision <- function(matrix) {
  # True positive
  tp <- matrix[2, 2]
  # false positive
  fp <- matrix[1, 2]
  return (tp / (tp + fp))
}
rec <- recall <- function(matrix) {
  # true positive
  tp <- matrix[2, 2]# false positive
  fn <- matrix[2, 1]
  return (tp / (tp + fn))
}
prec <- precision(table_mat)
rec <- recall(table_mat)
prec
rec
f1 <- 2 * ((prec * rec) / (prec + rec))
f1

#Make a ROC curve
ROCRpred <- prediction(predict, data$Type)
ROCRperf <- performance(ROCRpred, 'tpr', 'fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2, 1.7))
ROCRperf <- performance(ROCRpred, 'prec', 'rec')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2, 1.7))

#Cross validation
cost <- function(r, pi=0) mean(abs(r-pi)> 0.5)
cv.glm(data, glm_1, cost=cost)

