# Data analysis for Albano et al., submitted to Ecological Monographs
# Title: Adaptation to environmental variation in the native and introduced ranges of a cosmopolitan plant

# STEPS
# 1. Data set-up
# 2a. Experiment data analysis - Q1
# 2b. Experiment data analysis - Q2
# 2c. Experiment data analysis - Q3
# 2d. Experiment data analysis - Q4
# 2e. Experiment data analysis - Q5
# 3. Supplemental analysis - Climatic distance model

# Set the working directory
setwd("")

# Load all required packages
library(dplyr) # for data organization
library(ggplot2) # for plotting results and creating final figures
library(glmmTMB) # for running mixed effects models
library(car) # for obtaining model results (Anova)
library(emmeans) # for post-hoc tests on model results


#***************************** 
# 1. DATA SET-UP
#*****************************

# Models are labeled using a set of two or three numbers (format: #.#. or #.#.#), with each number corresponding to a category of organization

# The first number referring to the set of models, corresponding to the research question number in the manuscript
    # 1 = Q1
    # 2 = Q2
    # 3 = Q3
    # 4 = Q4
    # 5 = Q5

# The second number in the model name represents the response variable being tested
# "1st" refers to variables only collected for the 1st growing season
# "Full" refers to variables that are summed across the full length of the experiment (both growing seasons)
    # 1 = Survived1st (whether or not a plant survived the 1st growing season)
    # 2 = FloweredFull (whether or not a plant flowered during the full experiment)
    # 3 = SeededFull (whether or not a plant produced seeds during the full experiment)
    # 4 = FlowerHeadNumberFull (number of flower heads produced in the full experiment), also abbreviated as FH
    # 5 = SeedSetMassFull (mass of all seeds produced in the full experiment), also abbreviated as SS
    # 6 = MaxArea1st (the maximum plant surface area in the 1st growing season), also abbreviated as MA
    # 7 = GrowthRate1st (the growth rate of plants in the 1st growing season), also abbreviated as GR
    # 8 = Herbivory1st (the % of leaf area consumed in the 1st growing season), also abbreviated as HR

##In model set 4, the third number in the model name represented the bioclimatic variable in the model
    # 1 = Mean Annual Temperature
    # 5 = Max Temperature of Warmest Month
    # 6 = Min Temperature of Coldest Month
    # 12 = Annual Precipitation
    # 13 = Precipitation of Wettest Month
    # 14 = Precipitation of Driest Month

# Load and attach full dataset for the transcontinental common garden experiment
CG_Data = read.csv("Albano_et_al_2025_EM_Experiment_Data.csv", header = TRUE, sep = ",")
attach(CG_Data)

# convert variables to factors
CG_Data$Garden<-as.factor(CG_Data$Garden)
CG_Data$Position<-as.factor(CG_Data$Position)
CG_Data$Population<-as.factor(CG_Data$Population)
CG_Data$Rep<-as.factor(CG_Data$Rep)
CG_Data$ContOrigin<-as.factor(CG_Data$ContOrigin)
CG_Data$ContGarden<-as.factor(CG_Data$ContGarden)
CG_Data$Cyanotype<-as.factor(CG_Data$Cyanotype)
CG_Data$Cyanogenesis<-as.factor(CG_Data$Cyanogenesis)
CG_Data$Ac_ac<-as.factor(CG_Data$Ac_ac)
CG_Data$Li_li<-as.factor(CG_Data$Li_li)

# FlowerHeadNumberFull is zero-inflated so it must be analyzed in a two-step process
# First analyze FH data as a binomial (flowered or not), then remove all zeroes and analyze non-zero FH numbers
CG_FH_NoZero=subset(CG_Data, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# SeedSetMassFull is zero-inflated so it must be analyzed in a two-step process
# First analyze SS data as a binomial (produced seeds or not), then remove all zeroes and analyze non-zero SS masses
CG_SS_NoZero=subset(CG_Data, !SeedSetMassFull=="0") |> 
  droplevels()

# Many of the response variables are right skewed, so must be transformed into order to not violate model assumptions
# Perform appropriate transformations for all quantitative variables below
log_FH_Full<-log(CG_FH_NoZero$FlowerHeadNumberFull)
log_SS_Full<-log(CG_SS_NoZero$SeedSetMassFull)
log_MA_1st<-log(CG_Data$MaxArea1st)
sqrt_GR_1st<-sqrt(CG_Data$GrowthRate1st)
sqrt_HR_1st<-sqrt(CG_Data$Herbivory1st)

# Produce the three colour palettes needed for figure construction
Palette1 <- c("black")
Palette2 <- c("darkblue","green4", "darkorange2", "maroon3")
Palette3 <- c("darkgray","tan4")

#***************************** 
# 2. EXPERIMENT DATA ANALYSIS - a) MODELS TO ADDRESS Q1
#*****************************

# Step 1: Run glmmTMBs for each of the 8 response variables in the order listed above
    # Continent of Origin, Garden Location, and Weighted Climatic Distance included as predictors, plus all two-way interactions
    # Population included as a random effect
# Step 2: Run model summaries and Anova to obtain test statistics and p values
# Step 3: Run emtrends to obtain model-predicted slopes for each Garden
# Step 4: Build figure for each response variable against WeightedCD, grouped by Garden
# Step 5: For Flower Head Number and Seed Set Mass, also obtain slopes based on raw regression lines between WeightedCD and response

model1.1<-glmmTMB(Survived1st~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_Data, family = "binomial", na.action=na.omit)
summary(model1.1)
Anova(model1.1, p.adjust.method=TRUE)
model1.1em<-emtrends(model1.1, specs = pairwise~Garden, var = "WeightedCD")
print(model1.1em)

Fig1.1<-CG_Data |>
  filter(!is.na(Survived1st)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(Surv.mean = mean(Survived1st)) |>
  ggplot(aes(x = WeightedCD,  y = Surv.mean, group = Garden, color = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., 0), ymax = ..ymax.., fill = Garden)) +
  ggtitle("A") +
  ylab("Proportion Survived") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "") + 
  scale_fill_manual(values = Palette2, name = "") + 
  theme_classic() +
  coord_fixed(ratio = 7.2115384615) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig1.1 <- Fig1.1 + theme(legend.position = "none")
Fig1.1

model1.2<-glmmTMB(FloweredFull~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_Data, family="binomial", na.action=na.omit)
summary(model1.2)
Anova(model1.2, p.adjust.method=TRUE)
model1.2em<-emtrends(model1.2, specs = pairwise~Garden, var = "WeightedCD")
print(model1.2em)

Fig1.2<-CG_Data |>
  filter(!is.na(FloweredFull)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(Flowered.mean = mean(FloweredFull)) |>
  ggplot(aes(x = WeightedCD,  y = Flowered.mean, color = Garden, fill = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., 0), ymax = ..ymax.., fill = Garden)) +
  ggtitle("B") +
  ylab("Proportion Flowered") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "") + 
  scale_fill_manual(values = Palette2, name = "") + 
  theme_classic() +
  coord_fixed(ratio = 7.2115284615) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig1.2 <- Fig1.2 + theme(legend.position = "none")
Fig1.2

model1.3<-glmmTMB(SeededFull~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_Data, family="binomial", na.action=na.omit)
summary(model1.3)
Anova(model1.3, p.adjust.method=TRUE)
model1.3em<-emtrends(model1.3, specs = pairwise~Garden, var = "WeightedCD")
print(model1.3em)

Fig1.3<-CG_Data |>
  filter(!is.na(SeededFull)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(Seeded.mean = mean(SeededFull)) |>
  ggplot(aes(x = WeightedCD,  y = Seeded.mean, color = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., 0), ymax = ..ymax.., fill = Garden)) +
  ggtitle("C") +
  ylab("Proportion Producing Seeds") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "") + 
  scale_fill_manual(values = Palette2, name = "") + 
  theme_classic() +
  coord_fixed(ratio = 7.2115384615) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig1.3 <- Fig1.3 + theme(legend.position = "none")
Fig1.3

model1.4<-glmmTMB(log_FH_Full~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_FH_NoZero, na.action=na.omit)
summary(model1.4)
Anova(model1.4, p.adjust.method=TRUE)
model1.4em<-emtrends(model1.4, specs = pairwise~Garden, var = "WeightedCD")
print(model1.4em)

slopes1.4 <- CG_FH_NoZero |>
  filter(!is.na(FlowerHeadNumberFull)) |>
  group_by(Garden) |>
  summarize(
    slope = cov(WeightedCD, FlowerHeadNumberFull) / var(WeightedCD)) |>
  mutate(slope = format(round(slope, 5), nsmall = 5))
slopes1.4

Fig1.4<-CG_FH_NoZero |>
  filter(!is.na(FlowerHeadNumberFull)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(FH.mean = mean(FlowerHeadNumberFull)) |>
  ggplot(aes(x = WeightedCD,  y = FH.mean, group = Garden, color = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -20), ymax = ..ymax.., fill = Garden)) +  # clip lower CI at 0
  annotate("text", size=5, x=3, y=880, label="L*", hjust=0) +
  annotate("text", size=5, x=3, y=800, label="D*", hjust=0) +
  annotate("text", size=5, x=3, y=720, label="L×D*", hjust=0) +
  ggtitle("A") +
  ylab("Number of Flower Heads") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-20,1020), breaks = c(0,200,400,600,800,1000), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "") + 
  scale_fill_manual(values = Palette2, name = "") + 
  theme_classic() +
  coord_fixed(ratio = 0.0072115385) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12))
Fig1.4 <- Fig1.4 + theme(legend.position = "none")
Fig1.4

model1.5<-glmmTMB(log_SS_Full~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_SS_NoZero, na.action=na.omit)
summary(model1.5)
Anova(model1.5, p.adjust.method=TRUE)
model1.5em<-emtrends(model1.5, pairwise~Garden, var = "WeightedCD")
print(model1.5em)

slopes1.5 <- CG_SS_NoZero |>
  filter(!is.na(SeedSetMassFull)) |>
  group_by(Garden) |>
  summarize(
    slope = cov(WeightedCD, SeedSetMassFull) / var(WeightedCD)) |>
  mutate(slope = format(round(slope, 5), nsmall = 5))
slopes1.5

Fig1.5<-CG_SS_NoZero |>
  filter(!is.na(SeedSetMassFull)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(SS.mean = mean(SeedSetMassFull)) |>
  ggplot(aes(x = WeightedCD,  y = SS.mean, group = Garden, color = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.5), ymax = ..ymax.., fill = Garden)) +
  annotate("text", size=5, x=3, y=22, label="L*", hjust=0) +
  annotate("text", size=5, x=3, y=20, label="D*", hjust=0) +
  annotate("text", size=5, x=3, y=18, label="L×D*", hjust=0) +
  ggtitle("B") +
  ylab("Seed Set Mass (g)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.5,25.5), breaks = c(0,5,10,15,20,25), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "Garden Location", labels = c("Lafayette", "Mississauga", "Montpellier", "Uppsala")) + 
  scale_fill_manual(values = Palette2, name = "Garden Location", labels = c("Lafayette", "Mississauga", "Montpellier", "Uppsala")) + 
  theme_classic() +
  coord_fixed(ratio = 0.2884615385) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12))
Fig1.5

model1.6<-glmmTMB(log_MA_1st~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_Data, na.action=na.omit)
summary(model1.6)
Anova(model1.6, p.adjust.method=TRUE)
model1.6em<-emtrends(model1.6, specs = pairwise~Garden, var = "WeightedCD")
print(model1.6em)

Fig1.6<-CG_Data |>
  filter(!is.na(MaxArea1st)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(MA.mean = mean(MaxArea1st)) |>
  ggplot(aes(x = WeightedCD,  y = MA.mean, group = Garden, color = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -40), ymax = ..ymax.., fill = Garden)) +
  ggtitle("D") +
  ylab("Maximum Plant Area (cm²)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-40,2040), breaks = c(0,500,1000,1500,2000), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "Garden Location", labels = c("Lafayette", "Mississauga", "Montpellier", "Uppsala")) + 
  scale_fill_manual(values = Palette2, name = "Garden Location", labels = c("Lafayette", "Mississauga", "Montpellier", "Uppsala")) + 
  theme_classic() +
  coord_fixed(ratio = 0.0036057692) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.position=c(0.65,0.75))
Fig1.6

model1.7<-glmmTMB(sqrt_GR_1st~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_Data, na.action=na.omit)
summary(model1.7)
Anova(model1.7, p.adjust.method=TRUE)
model1.7em<-emtrends(model1.7, specs = pairwise~Garden, var = "WeightedCD")
print(model1.7em)

Fig1.7<-CG_Data |>
  filter(!is.na(GrowthRate1st)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(GR.mean = mean(GrowthRate1st)) |>
  ggplot(aes(x = WeightedCD,  y = GR.mean, group = Garden, color = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., 0), ymax = ..ymax.., fill = Garden)) +
  ggtitle("E") +
  ylab("Growth Rate (cm²/day)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.003,0.153), breaks = c(0,0.05,0.10,0.15), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "") + 
  scale_fill_manual(values = Palette2, name = "") + 
  theme_classic() +
  coord_fixed(ratio = 48.0769230769) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig1.7 <- Fig1.7 + theme(legend.position = "none")
Fig1.7

model1.8<-glmmTMB(sqrt_HR_1st~ContOrigin+Garden+WeightedCD+ContOrigin:Garden+ContOrigin:WeightedCD+Garden:WeightedCD+(1|Population), data=CG_Data, na.action=na.omit)
summary(model1.8)
Anova(model1.8, p.adjust.method=TRUE)
model1.8em<-emtrends(model1.8, specs = pairwise~Garden, var = "WeightedCD")
print(model1.8em)

Fig1.8<-CG_Data |>
  filter(!is.na(Herbivory1st)) |>
  group_by(WeightedCD, Garden) |>
  dplyr::summarize(HR.mean = mean(Herbivory1st)) |>
  ggplot(aes(x = WeightedCD,  y = HR.mean, group = Garden, color = Garden)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., 0), ymax = ..ymax.., fill = Garden)) +
  ggtitle("F") +
  ylab("% Leaf Area Consumed") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-1.5,76.5), breaks = c(0,25,50,75), expand = c(0,0)) +
  scale_color_manual(values = Palette2, name = "") + 
  scale_fill_manual(values = Palette2, name = "") + 
  theme_classic() +
  coord_fixed(ratio = 0.0961538462) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12))
Fig1.8 <- Fig1.8 + theme(legend.position = "none")
Fig1.8

#***************************** 
# 2. EXPERIMENT DATA ANALYSIS - b) MODELS TO ADDRESS Q2
#*****************************

# Step 1: Run glmmTMBs for each of the 8 response variables in the order listed above
# Continent of Origin, Continent of Garden, and Weighted Climatic Distance included as predictors, plus all two- and three-way interactions
# Population included as a random effect
# Step 2: Run model summaries and Anova to obtain test statistics and p values
# Step 3: Run emtrends to obtain model-predicted slopes for each ContOrigin/ContGarden combination
# Step 4: Build figure for each response variable against WeightedCD, grouped by Continent of Origin and plotted separately by Continent of Garden
# Step 5: For Flower Head Number and Seed Set Mass, also obtain slopes based on raw regression lines between WeightedCD and response

model2.1<-glmmTMB(Survived1st~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_Data, family = "binomial", na.action=na.omit)
summary(model2.1)
Anova(model2.1, p.adjust.method=TRUE)
model2.1em<-emtrends(model2.1, specs = pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.1em)

Fig2.1.NA<-CG_Data |>
  filter(ContGarden == "America") |>
  filter(!is.na(Survived1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(Surv.mean = mean(Survived1st)) |>
  ggplot(aes(x = WeightedCD,  y = Surv.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("A") +
  ylab("Proportion Survived") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 9.6153846154) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.1.NA <- Fig2.1.NA + theme(legend.position = "none")
Fig2.1.NA

Fig2.1.EU<-CG_Data |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(Survived1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(Surv.mean = mean(Survived1st)) |>
  ggplot(aes(x = WeightedCD,  y = Surv.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("B") +
  ylab("Proportion Survived") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 9.6153846154) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.1.EU <- Fig2.1.EU + theme(legend.position = "none")
Fig2.1.EU

model2.2<-glmmTMB(FloweredFull~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_Data, family = "binomial", na.action=na.omit)
summary(model2.2)
Anova(model2.2, p.adjust.method=TRUE)
model2.2em<-emtrends(model2.2, specs = pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.2em)

Fig2.2.NA<-CG_Data |>
  filter(ContGarden == "America") |>
  filter(!is.na(FloweredFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(Flowered.mean = mean(FloweredFull)) |>
  ggplot(aes(x = WeightedCD,  y = Flowered.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("C") +
  ylab("Proportion Flowered") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 9.6153846154) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.2.NA <- Fig2.2.NA + theme(legend.position = "none")
Fig2.2.NA

Fig2.2.EU<-CG_Data |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(FloweredFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(Flowered.mean = mean(FloweredFull)) |>
  ggplot(aes(x = WeightedCD,  y = Flowered.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("D") +
  ylab("Proportion Flowered") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 9.6153846154) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.2.EU <- Fig2.2.EU + theme(legend.position = "none")
Fig2.2.EU

model2.3<-glmmTMB(SeededFull~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_Data, family = "binomial", na.action=na.omit)
summary(model2.3)
Anova(model2.3, p.adjust.method=TRUE)
model2.3em<-emtrends(model2.3, specs = pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.3em)

Fig2.3.NA<-CG_Data |>
  filter(ContGarden == "America") |>
  filter(!is.na(SeededFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(Seeded.mean = mean(SeededFull)) |>
  ggplot(aes(x = WeightedCD,  y = Seeded.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("E") +
  ylab("Proportion Producing Seeds") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 9.6153846154) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.3.NA <- Fig2.3.NA + theme(legend.position = "none")
Fig2.3.NA

Fig2.3.EU<-CG_Data |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(SeededFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(Seeded.mean = mean(SeededFull)) |>
  ggplot(aes(x = WeightedCD,  y = Seeded.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("F") +
  ylab("Proportion Producing Seeds") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 9.6153846154) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.3.EU <- Fig2.3.EU + theme(legend.position = "none")
Fig2.3.EU

model2.4<-glmmTMB(log_FH_Full~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_FH_NoZero, na.action=na.omit)
summary(model2.4)
Anova(model2.4, p.adjust.method=TRUE)
model2.4em<-emtrends(model2.4, pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.4em)

slopes2.4 <- CG_FH_NoZero |>
  filter(!is.na(FlowerHeadNumberFull)) |>
  group_by(ContGarden, ContOrigin) |>
  summarize(
    slope = cov(WeightedCD, FlowerHeadNumberFull) / var(WeightedCD)) |>
  mutate(slope = format(round(slope, 5), nsmall = 5))
slopes2.4

Fig2.4.NA<-CG_FH_NoZero |>
  filter(ContGarden == "America") |>
  filter(!is.na(FlowerHeadNumberFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(FH.mean = mean(FlowerHeadNumberFull)) |>
  ggplot(aes(x = WeightedCD,  y = FH.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -12), ymax = ..ymax.., fill = ContOrigin)) +
  annotate("text", size=5, x=3, y=528, label="D***", hjust=0) +
  annotate("text", size=5, x=3, y=480, label="O×G×D***", hjust=0) +
  ggtitle("             North American Gardens
A") +
  ylab("Number of Flower Heads") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-12,612), breaks = c(0,200,400,600), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 0.0120192308) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.4.NA <- Fig2.4.NA + theme(legend.position = "none")
Fig2.4.NA

Fig2.4.EU<-CG_FH_NoZero |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(FlowerHeadNumberFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(FH.mean = mean(FlowerHeadNumberFull)) |>
  ggplot(aes(x = WeightedCD,  y = FH.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -20), ymax = ..ymax.., fill = ContOrigin)) +
  annotate("text", size=5, x=3, y=880, label="D***", hjust=0) +
  annotate("text", size=5, x=3, y=800, label="O×G×D**", hjust=0) +
  ggtitle("               European Gardens
B") +
  ylab("Number of Flower Heads") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-20,1020), breaks = c(0,200,400,600,800,1000), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 0.0072115385) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.4.EU <- Fig2.4.EU + theme(legend.position = "none")
Fig2.4.EU

model2.5<-glmmTMB(log_SS_Full~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_SS_NoZero, na.action=na.omit)
summary(model2.5)
Anova(model2.5, p.adjust.method=TRUE)
model2.5em<-emtrends(model2.5, pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.5em)

slopes2.5 <- CG_SS_NoZero |>
  filter(!is.na(SeedSetMassFull)) |>
  group_by(ContGarden, ContOrigin) |>
  summarize(
    slope = cov(WeightedCD, SeedSetMassFull) / var(WeightedCD)) |>
  mutate(slope = format(round(slope, 5), nsmall = 5))
slopes2.5

Fig2.5.NA<-CG_SS_NoZero |>
  filter(ContGarden == "America") |>
  filter(!is.na(SeedSetMassFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(SS.mean = mean(SeedSetMassFull)) |>
  ggplot(aes(x = WeightedCD,  y = SS.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.4), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("
C") +
  ylab("Seed Set Mass (g)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.4,20.4), breaks = c(0,5,10,15,20), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 0.3605769231) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.5.NA <- Fig2.5.NA + theme(legend.position = "none")
Fig2.5.NA

Fig2.5.EU<-CG_SS_NoZero |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(SeedSetMassFull)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(SS.mean = mean(SeedSetMassFull)) |>
  ggplot(aes(x = WeightedCD,  y = SS.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.5), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("
D") +
  ylab("Seed Set Mass (g)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.5,25.5), breaks = c(0,5,10,15,20,25), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "Continent of Origin", labels = c("Europe", "North America")) + 
  scale_fill_manual(values = Palette3, name = "Continent of Origin", labels = c("Europe", "North America")) + 
  coord_fixed(ratio = 0.2884615385) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.position=c(0.65,0.65))
Fig2.5.EU

model2.6<-glmmTMB(log_MA_1st~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_Data, na.action=na.omit)
summary(model2.6)
Anova(model2.6, p.adjust.method=TRUE)
model2.6em<-emtrends(model2.6, specs = pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.6em)

Fig2.6.NA<-CG_Data |>
  filter(ContGarden == "America") |>
  filter(!is.na(MaxArea1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(MA.mean = mean(MaxArea1st)) |>
  ggplot(aes(x = WeightedCD,  y = MA.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -16), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("G") +
  ylab("Maximum Plant Area (cm²)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-16,816), breaks = c(0,200,400,600,800), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 0.0120192308) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.6.NA <- Fig2.6.NA + theme(legend.position = "none")
Fig2.6.NA

Fig2.6.EU<-CG_Data |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(MaxArea1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(MA.mean = mean(MaxArea1st)) |>
  ggplot(aes(x = WeightedCD,  y = MA.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -40), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("H") +
  ylab("Maximum Plant Area (cm²)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-40,2040), breaks = c(0,500,1000,1500,2000), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 0.0048076923) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.6.EU <- Fig2.6.EU + theme(legend.position = "none")
Fig2.6.EU

model2.7<-glmmTMB(sqrt_GR_1st~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_Data, na.action=na.omit)
summary(model2.7)
Anova(model2.7, p.adjust.method=TRUE)
model2.7em<-emtrends(model2.7, specs = pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.7em)

Fig2.7.NA<-CG_Data |>
  filter(ContGarden == "America") |>
  filter(!is.na(GrowthRate1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(GR.mean = mean(GrowthRate1st)) |>
  ggplot(aes(x = WeightedCD,  y = GR.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., 0), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("I") +
  ylab("Growth Rate (cm²/day)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.003,0.153), breaks = c(0,0.05,0.10,0.15), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") +
  coord_fixed(ratio = 64.1025641026) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.7.NA <- Fig2.7.NA + theme(legend.position = "none")
Fig2.7.NA

Fig2.7.EU<-CG_Data |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(GrowthRate1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(GR.mean = mean(GrowthRate1st)) |>
  ggplot(aes(x = WeightedCD,  y = GR.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., 0), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("J") +
  ylab("Growth Rate (cm²/day)") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-0.002,0.102), breaks = c(0,0.05,0.10), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 96.1538461538) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.7.EU <- Fig2.7.EU + theme(legend.position = "none")
Fig2.7.EU

model2.8<-glmmTMB(sqrt_HR_1st~ContOrigin*ContGarden*WeightedCD+(1|Population), data=CG_Data, na.action=na.omit)
summary(model2.8)
Anova(model2.8)
model2.8em<-emtrends(model2.8, specs = pairwise~ContOrigin*ContGarden, var = "WeightedCD")
print(model2.8em)

Fig2.8.NA<-CG_Data |>
  filter(ContGarden == "America") |>
  filter(!is.na(Herbivory1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(HR.mean = mean(Herbivory1st)) |>
  ggplot(aes(x = WeightedCD,  y = HR.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -1.5), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("K") +
  ylab("% Leaf Area Consumed") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-1.5,76.5), breaks = c(0,25,50,75,100), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 0.1282051282) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.8.NA <- Fig2.8.NA + theme(legend.position = "none")
Fig2.8.NA

Fig2.8.EU<-CG_Data |>
  filter(ContGarden == "Europe") |>
  filter(!is.na(Herbivory1st)) |>
  group_by(WeightedCD, ContOrigin) |>
  dplyr::summarize(HR.mean = mean(Herbivory1st)) |>
  ggplot(aes(x = WeightedCD,  y = HR.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -1.5), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("L") +
  ylab("% Leaf Area Consumed") +
  xlab("Climatic Distance") +
  scale_y_continuous(limits = c(-1.5,76.5), breaks = c(0,25,50,75,100), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "") + 
  scale_fill_manual(values = Palette3, name = "") + 
  coord_fixed(ratio = 0.1282051282) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig2.8.EU <- Fig2.8.EU + theme(legend.position = "none")
Fig2.8.EU

#***************************** 
# EXPERIMENT DATA ANALYSIS - c) MODELS TO ADDRESS Q3
#*****************************

# Step 1: Analyze and plot Ac/ac and Li/li data for presence of latitudinal clines based on ContOrigin

# Subset data for the full experiment down to only the individuals with cyanogenesis data present (remove NAs)
CG_Cyano=subset(CG_Data, !Cyanotype=="NA") |> 
  droplevels()

##Use glmmTMBs to test for the presence of clines in Ac across the entire experiment
model3.Ac<-glmmTMB(Ac_ac_Num~Latitude*ContOrigin+(1|ContOrigin/Population), data=CG_Cyano, family = "binomial", na.action=na.omit)
summary(model3.Ac)
Anova(model3.Ac, p.adjust.method=TRUE)

Fig3.Ac<-CG_Cyano |>
  group_by(Latitude, ContOrigin) |>
  dplyr::summarize(Ac.mean = mean(Ac_ac_Num)) |>
  ggplot(aes(x = Latitude,  y = Ac.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("A") +
  ylab("Proportion Ac") +
  xlab("Latitude") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "Continent of Origin", labels = c("Europe", "North America")) + 
  scale_fill_manual(values = Palette3, name = "Continent of Origin", labels = c("Europe", "North America")) + 
  theme_classic() +
  coord_fixed(ratio = 40) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12))
Fig3.Ac <- Fig3.Ac + theme(legend.position = "none")
Fig3.Ac

##Use glmmTMBs to test for the presence of clines in Li across the entire experiment
model3.Li<-glmmTMB(Li_li_Num~Latitude*ContOrigin+(1|ContOrigin/Population), data=CG_Cyano, family = "binomial", na.action=na.omit)
summary(model3.Li)
Anova(model3.Li, p.adjust.method=TRUE)

Fig3.Li<-CG_Cyano |>
  group_by(Latitude, ContOrigin) |>
  dplyr::summarize(Li.mean = mean(Li_li_Num)) |>
  ggplot(aes(x = Latitude,  y = Li.mean, group = ContOrigin, color = ContOrigin)) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax.., fill = ContOrigin)) +
  ggtitle("B") +
  ylab("Proportion Li") +
  xlab("Latitude") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette3, name = "Continent of Origin", labels = c("Europe", "North America")) + 
  scale_fill_manual(values = Palette3, name = "Continent of Origin", labels = c("Europe", "North America")) + 
  theme_classic() +
  coord_fixed(ratio = 40) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12))
Fig3.Li

# Step 2: Set up fitness and herbivory data for analysis

# FlowerHeadNumberFull is zero-inflated so it must be analyzed in a two-step process
# First analyze FH data as a binomial (flowered or not), then remove all zeroes and analyze non-zero FH numbers
CG_Cyano_FH_NoZero=subset(CG_Cyano, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# SeedSetMassFull is zero-inflated so it must be analyzed in a two-step process
# First analyze SS data as a binomial (produced seeds or not), then remove all zeroes and analyze non-zero SS numbers
CG_Cyano_SS_NoZero=subset(CG_Cyano, !SeedSetMassFull=="0") |> 
  droplevels()

# Many of the response variables are right skewed, so must be transformed into order to not violate model assumptions
# Perform appropriate transformations for all quantitative variables below
log_FH_Cyano<-log(CG_Cyano_FH_NoZero$FlowerHeadNumberFull)
log_SS_Cyano<-log(CG_Cyano_SS_NoZero$SeedSetMassFull)
log_MA_Cyano<-log(CG_Cyano$MaxArea1st)
log_GR_Cyano<-log(CG_Cyano$GrowthRate1st)
sqrt_HR_Cyano<-sqrt(CG_Cyano$Herbivory1st)

# Step 3: Run glmmTMBs for each of the 8 response variables in the order listed above
    # Cyanotype, Garden Location, and Weighted Climatic Distance included as predictors, plus all interactions that include Cyanotype
    # Population included as a random effect
# Step 4: Run model summaries and Anova to obtain test statistics and p values

model3.1<-glmmTMB(Survived1st~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano, family = "binomial", na.action=na.omit)
summary(model3.1)
Anova(model3.1, p.adjust.method=TRUE)

model3.2<-glmmTMB(FloweredFull~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano, family = "binomial", na.action=na.omit)
summary(model3.2)
Anova(model3.2, p.adjust.method=TRUE)

model3.3<-glmmTMB(SeededFull~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano, family = "binomial", na.action=na.omit)
summary(model3.3)
Anova(model3.3, p.adjust.method=TRUE)

model3.4<-glmmTMB(log_FH_Cyano~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano_FH_NoZero, na.action=na.omit)
summary(model3.4)
Anova(model3.4, p.adjust.method=TRUE)

model3.5<-glmmTMB(log_SS_Cyano~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano_SS_NoZero, na.action=na.omit)
summary(model3.5)
Anova(model3.5, p.adjust.method=TRUE)

model3.6<-glmmTMB(log_MA_Cyano~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano, na.action=na.omit)
summary(model3.6)
Anova(model3.6, p.adjust.method=TRUE)

model3.7<-glmmTMB(log_GR_Cyano~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano, na.action=na.omit)
summary(model3.7)
Anova(model3.7, p.adjust.method=TRUE)

model3.8<-glmmTMB(sqrt_HR_Cyano~Cyanotype+Garden+WeightedCD+
                    Cyanotype:Garden+Cyanotype:WeightedCD+Cyanotype:Garden:WeightedCD+
                    (1|Population), data=CG_Cyano, na.action=na.omit)
summary(model3.8)
Anova(model3.8, p.adjust.method=TRUE)

#***************************** 
# 2. EXPERIMENT DATA ANALYSIS - d) MODELS TO ADDRESS Q4
#*****************************

# Step 1: Run glmmTMBs for each of the 7 response variables in the order listed above (excluding Herbivory)
    # Each of the 6 bioclimatic variables (both linear and quadratic terms incuded as predictors) run in a separate model with the number system listed:
      # 1 = Mean Annual Temperature
      # 5 = Max Temperature of Warmest Month
      # 6 = Min Temperature of Coldest Month
      # 12 = Annual Precipitation
      # 13 = Precipitation of Wettest Month
      # 14 = Precipitation of Driest Month
    # Population included as a random effect
# Step 2: Run model summaries and Anova to obtain test statistics and p values

#***************************** 
# Q4a. Mississauga
#*****************************

# Subset full experiment data to just data from the Mississauga common garden
CG_Mis=subset(CG_Data, Garden=="Mississauga")

# Subset CG_Mis to analyze zero-inflated flower head numbers as above
CG_Mis_FH_NoZero=subset(CG_Mis, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# Subset CG_Mis to analyze zero-inflated seed set masses as above
CG_Mis_SS_NoZero=subset(CG_Mis, !SeedSetMassFull=="0") |> 
  droplevels()

# Perform appropriate transformations for all quantitative variables below
# Each quantitative variable is also used in its z-score standardized form for equal comparisons across all models
log_FH_Mis<-log(CG_Mis_FH_NoZero$FH_Full_z)
log_SS_Mis<-log(CG_Mis_SS_NoZero$SS_Full_z)
log_MA_Mis<-log(CG_Mis$MA_1st_z)

# Q4a - BIO1
model4.1.1_Mis<-glmmTMB(Survived1st~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.1.1_Mis)
Anova(model4.1.1_Mis, p.adjust.method=TRUE)

model4.2.1_Mis<-glmmTMB(FloweredFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.2.1_Mis)
Anova(model4.2.1_Mis, p.adjust.method=TRUE)

model4.3.1_Mis<-glmmTMB(SeededFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.3.1_Mis)
Anova(model4.3.1_Mis, p.adjust.method=TRUE)

model4.4.1_Mis<-glmmTMB(log_FH_Mis~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mis_FH_NoZero, na.action=na.omit)
summary(model4.4.1_Mis)
Anova(model4.4.1_Mis, p.adjust.method=TRUE)

model4.5.1_Mis<-glmmTMB(log_SS_Mis~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mis_SS_NoZero, na.action=na.omit)
summary(model4.5.1_Mis)
Anova(model4.5.1_Mis, p.adjust.method=TRUE)

model4.6.1_Mis<-glmmTMB(log_MA_Mis~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.6.1_Mis)
Anova(model4.6.1_Mis, p.adjust.method=TRUE)

model4.7.1_Mis<-glmmTMB(GR_1st_z~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.7.1_Mis)
Anova(model4.7.1_Mis, p.adjust.method=TRUE)

# Q4a - BIO5
model4.1.5_Mis<-glmmTMB(Survived1st~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.1.5_Mis)
Anova(model4.1.5_Mis, p.adjust.method=TRUE)

model4.2.5_Mis<-glmmTMB(FloweredFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.2.5_Mis)
Anova(model4.2.5_Mis, p.adjust.method=TRUE)

model4.3.5_Mis<-glmmTMB(SeededFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.3.5_Mis)
Anova(model4.3.5_Mis, p.adjust.method=TRUE)

model4.4.5_Mis<-glmmTMB(log_FH_Mis~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mis_FH_NoZero, na.action=na.omit)
summary(model4.4.5_Mis)
Anova(model4.4.5_Mis, p.adjust.method=TRUE)

model4.5.5_Mis<-glmmTMB(log_SS_Mis~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mis_SS_NoZero, na.action=na.omit)
summary(model4.5.5_Mis)
Anova(model4.5.5_Mis, p.adjust.method=TRUE)

model4.6.5_Mis<-glmmTMB(log_MA_Mis~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.6.5_Mis)
Anova(model4.6.5_Mis, p.adjust.method=TRUE)

model4.7.5_Mis<-glmmTMB(GR_1st_z~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.7.5_Mis)
Anova(model4.7.5_Mis, p.adjust.method=TRUE)

# Q4a - BIO6
model4.1.6_Mis<-glmmTMB(Survived1st~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.1.6_Mis)
Anova(model4.1.6_Mis, p.adjust.method=TRUE)

model4.2.6_Mis<-glmmTMB(FloweredFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.2.6_Mis)
Anova(model4.2.6_Mis, p.adjust.method=TRUE)

model4.3.6_Mis<-glmmTMB(SeededFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.3.6_Mis)
Anova(model4.3.6_Mis, p.adjust.method=TRUE)

model4.4.6_Mis<-glmmTMB(log_FH_Mis~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mis_FH_NoZero, na.action=na.omit)
summary(model4.4.6_Mis)
Anova(model4.4.6_Mis, p.adjust.method=TRUE)

model4.5.6_Mis<-glmmTMB(log_SS_Mis~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mis_SS_NoZero, na.action=na.omit)
summary(model4.5.6_Mis)
Anova(model4.5.6_Mis, p.adjust.method=TRUE)

model4.6.6_Mis<-glmmTMB(log_MA_Mis~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.6.6_Mis)
Anova(model4.6.6_Mis, p.adjust.method=TRUE)

model4.7.6_Mis<-glmmTMB(GR_1st_z~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.7.6_Mis)
Anova(model4.7.6_Mis, p.adjust.method=TRUE)

# Q4a - BIO12
model4.1.12_Mis<-glmmTMB(Survived1st~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.1.12_Mis)
Anova(model4.1.12_Mis, p.adjust.method=TRUE)

model4.2.12_Mis<-glmmTMB(FloweredFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.2.12_Mis)
Anova(model4.2.12_Mis, p.adjust.method=TRUE)

model4.3.12_Mis<-glmmTMB(SeededFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.3.12_Mis)
Anova(model4.3.12_Mis, p.adjust.method=TRUE)

model4.4.12_Mis<-glmmTMB(log_FH_Mis~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mis_FH_NoZero, na.action=na.omit)
summary(model4.4.12_Mis)
Anova(model4.4.12_Mis, p.adjust.method=TRUE)

model4.5.12_Mis<-glmmTMB(log_SS_Mis~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mis_SS_NoZero, na.action=na.omit)
summary(model4.5.12_Mis)
Anova(model4.5.12_Mis, p.adjust.method=TRUE)

model4.6.12_Mis<-glmmTMB(log_MA_Mis~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.6.12_Mis)
Anova(model4.6.12_Mis, p.adjust.method=TRUE)

model4.7.12_Mis<-glmmTMB(GR_1st_z~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.7.12_Mis)
Anova(model4.7.12_Mis, p.adjust.method=TRUE)

# Q4a - BIO13
model4.1.13_Mis<-glmmTMB(Survived1st~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.1.13_Mis)
Anova(model4.1.13_Mis, p.adjust.method=TRUE)

model4.2.13_Mis<-glmmTMB(FloweredFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.2.13_Mis)
Anova(model4.2.13_Mis, p.adjust.method=TRUE)

model4.3.13_Mis<-glmmTMB(SeededFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.3.13_Mis)
Anova(model4.3.13_Mis, p.adjust.method=TRUE)

model4.4.13_Mis<-glmmTMB(log_FH_Mis~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mis_FH_NoZero, na.action=na.omit)
summary(model4.4.13_Mis)
Anova(model4.4.13_Mis, p.adjust.method=TRUE)

model4.5.13_Mis<-glmmTMB(log_SS_Mis~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mis_SS_NoZero, na.action=na.omit)
summary(model4.5.13_Mis)
Anova(model4.5.13_Mis, p.adjust.method=TRUE)

model4.6.13_Mis<-glmmTMB(log_MA_Mis~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.6.13_Mis)
Anova(model4.6.13_Mis, p.adjust.method=TRUE)

model4.7.13_Mis<-glmmTMB(GR_1st_z~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.7.13_Mis)
Anova(model4.7.13_Mis, p.adjust.method=TRUE)

# Q4a - BIO14
model4.1.14_Mis<-glmmTMB(Survived1st~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.1.14_Mis)
Anova(model4.1.14_Mis, p.adjust.method=TRUE)

model4.2.14_Mis<-glmmTMB(FloweredFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.2.14_Mis)
Anova(model4.2.14_Mis, p.adjust.method=TRUE)

model4.3.14_Mis<-glmmTMB(SeededFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mis, family = "binomial", na.action=na.omit)
summary(model4.3.14_Mis)
Anova(model4.3.14_Mis, p.adjust.method=TRUE)

model4.4.14_Mis<-glmmTMB(log_FH_Mis~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mis_FH_NoZero, na.action=na.omit)
summary(model4.4.14_Mis)
Anova(model4.4.14_Mis, p.adjust.method=TRUE)

model4.5.14_Mis<-glmmTMB(log_SS_Mis~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mis_SS_NoZero, na.action=na.omit)
summary(model4.5.14_Mis)
Anova(model4.5.14_Mis, p.adjust.method=TRUE)

model4.6.14_Mis<-glmmTMB(log_MA_Mis~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.6.14_Mis)
Anova(model4.6.14_Mis, p.adjust.method=TRUE)

model4.7.14_Mis<-glmmTMB(GR_1st_z~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mis, na.action=na.omit)
summary(model4.7.14_Mis)
Anova(model4.7.14_Mis, p.adjust.method=TRUE)

#***************************** 
# Q4b. Lafayette
#*****************************

# Subset full experiment data to just data from the Lafayette common garden
CG_Laf=subset(CG_Data, Garden=="Lafayette")

# Subset CG_Laf to analyze zero-inflated flower head numbers as above
CG_Laf_FH_NoZero=subset(CG_Laf, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# Subset CG_Laf to analyze zero-inflated seed set masses as above
CG_Laf_SS_NoZero=subset(CG_Laf, !SeedSetMassFull=="0") |> 
  droplevels()

# Perform appropriate transformations for all quantitative variables below
# Each quantitative variable is also used in its z-score standardized form for equal comparisons across all models
log_FH_Laf<-log(CG_Laf_FH_NoZero$FH_Full_z)
sqrt_SS_Laf<-sqrt(CG_Laf_SS_NoZero$SS_Full_z)
log_MA_Laf<-log(CG_Laf$MA_1st_z)
log_GR_Laf<-log(CG_Laf$GR_1st_z)

# Q4b. BIO1
model4.1.1_Laf<-glmmTMB(Survived1st~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.1.1_Laf)
Anova(model4.1.1_Laf, p.adjust.method=TRUE)

model4.2.1_Laf<-glmmTMB(FloweredFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.2.1_Laf)
Anova(model4.2.1_Laf, p.adjust.method=TRUE)

model4.3.1_Laf<-glmmTMB(SeededFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.3.1_Laf)
Anova(model4.3.1_Laf, p.adjust.method=TRUE)

model4.4.1_Laf<-glmmTMB(log_FH_Laf~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Laf_FH_NoZero, na.action=na.omit)
summary(model4.4.1_Laf)
Anova(model4.4.1_Laf, p.adjust.method=TRUE)

model4.5.1_Laf<-glmmTMB(sqrt_SS_Laf~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Laf_SS_NoZero, na.action=na.omit)
summary(model4.5.1_Laf)
Anova(model4.5.1_Laf, p.adjust.method=TRUE)

model4.6.1_Laf<-glmmTMB(log_MA_Laf~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.6.1_Laf)
Anova(model4.6.1_Laf, p.adjust.method=TRUE)

model4.7.1_Laf<-glmmTMB(log_GR_Laf~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.7.1_Laf)
Anova(model4.7.1_Laf, p.adjust.method=TRUE)

# Q4b. BIO5
model4.1.5_Laf<-glmmTMB(Survived1st~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.1.5_Laf)
Anova(model4.1.5_Laf, p.adjust.method=TRUE)

model4.2.5_Laf<-glmmTMB(FloweredFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.2.5_Laf)
Anova(model4.2.5_Laf, p.adjust.method=TRUE)

model4.3.5_Laf<-glmmTMB(SeededFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.3.5_Laf)
Anova(model4.3.5_Laf, p.adjust.method=TRUE)

model4.4.5_Laf<-glmmTMB(log_FH_Laf~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Laf_FH_NoZero, na.action=na.omit)
summary(model4.4.5_Laf)
Anova(model4.4.5_Laf, p.adjust.method=TRUE)

model4.5.5_Laf<-glmmTMB(sqrt_SS_Laf~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Laf_SS_NoZero, na.action=na.omit)
summary(model4.5.5_Laf)
Anova(model4.5.5_Laf, p.adjust.method=TRUE)

model4.6.5_Laf<-glmmTMB(log_MA_Laf~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.6.5_Laf)
Anova(model4.6.5_Laf, p.adjust.method=TRUE)

model4.7.5_Laf<-glmmTMB(log_GR_Laf~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.7.5_Laf)
Anova(model4.7.5_Laf, p.adjust.method=TRUE)

# Q4b. BIO6
model4.1.6_Laf<-glmmTMB(Survived1st~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.1.6_Laf)
Anova(model4.1.6_Laf, p.adjust.method=TRUE)

model4.2.6_Laf<-glmmTMB(FloweredFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.2.6_Laf)
Anova(model4.2.6_Laf, p.adjust.method=TRUE)

model4.3.6_Laf<-glmmTMB(SeededFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.3.6_Laf)
Anova(model4.3.6_Laf, p.adjust.method=TRUE)

model4.4.6_Laf<-glmmTMB(log_FH_Laf~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Laf_FH_NoZero, na.action=na.omit)
summary(model4.4.6_Laf)
Anova(model4.4.6_Laf, p.adjust.method=TRUE)

model4.5.6_Laf<-glmmTMB(sqrt_SS_Laf~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Laf_SS_NoZero, na.action=na.omit)
summary(model4.5.6_Laf)
Anova(model4.5.6_Laf, p.adjust.method=TRUE)

model4.6.6_Laf<-glmmTMB(log_MA_Laf~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.6.6_Laf)
Anova(model4.6.6_Laf, p.adjust.method=TRUE)

model4.7.6_Laf<-glmmTMB(log_GR_Laf~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.7.6_Laf)
Anova(model4.7.6_Laf, p.adjust.method=TRUE)

# Q4b. BIO12
model4.1.12_Laf<-glmmTMB(Survived1st~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.1.12_Laf)
Anova(model4.1.12_Laf, p.adjust.method=TRUE)

model4.2.12_Laf<-glmmTMB(FloweredFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.2.12_Laf)
Anova(model4.2.12_Laf, p.adjust.method=TRUE)

model4.3.12_Laf<-glmmTMB(SeededFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.3.12_Laf)
Anova(model4.3.12_Laf, p.adjust.method=TRUE)

model4.4.12_Laf<-glmmTMB(log_FH_Laf~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Laf_FH_NoZero, na.action=na.omit)
summary(model4.4.12_Laf)
Anova(model4.4.12_Laf, p.adjust.method=TRUE)

model4.5.12_Laf<-glmmTMB(sqrt_SS_Laf~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Laf_SS_NoZero, na.action=na.omit)
summary(model4.5.12_Laf)
Anova(model4.5.12_Laf, p.adjust.method=TRUE)

model4.6.12_Laf<-glmmTMB(log_MA_Laf~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.6.12_Laf)
Anova(model4.6.12_Laf, p.adjust.method=TRUE)

model4.7.12_Laf<-glmmTMB(log_GR_Laf~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.7.12_Laf)
Anova(model4.7.12_Laf, p.adjust.method=TRUE)

# Q4b. BIO13
model4.1.13_Laf<-glmmTMB(Survived1st~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.1.13_Laf)
Anova(model4.1.13_Laf, p.adjust.method=TRUE)

model4.2.13_Laf<-glmmTMB(FloweredFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.2.13_Laf)
Anova(model4.2.13_Laf, p.adjust.method=TRUE)

model4.3.13_Laf<-glmmTMB(SeededFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.3.13_Laf)
Anova(model4.3.13_Laf, p.adjust.method=TRUE)

model4.4.13_Laf<-glmmTMB(log_FH_Laf~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Laf_FH_NoZero, na.action=na.omit)
summary(model4.4.13_Laf)
Anova(model4.4.13_Laf, p.adjust.method=TRUE)

model4.5.13_Laf<-glmmTMB(sqrt_SS_Laf~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Laf_SS_NoZero, na.action=na.omit)
summary(model4.5.13_Laf)
Anova(model4.5.13_Laf, p.adjust.method=TRUE)

model4.6.13_Laf<-glmmTMB(log_MA_Laf~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.6.13_Laf)
Anova(model4.6.13_Laf, p.adjust.method=TRUE)

model4.7.13_Laf<-glmmTMB(log_GR_Laf~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.7.13_Laf)
Anova(model4.7.13_Laf, p.adjust.method=TRUE)

# Q4b. BIO14
model4.1.14_Laf<-glmmTMB(Survived1st~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.1.14_Laf)
Anova(model4.1.14_Laf, p.adjust.method=TRUE)

model4.2.14_Laf<-glmmTMB(FloweredFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.2.14_Laf)
Anova(model4.2.14_Laf, p.adjust.method=TRUE)

model4.3.14_Laf<-glmmTMB(SeededFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Laf, family = "binomial", na.action=na.omit)
summary(model4.3.14_Laf)
Anova(model4.3.14_Laf, p.adjust.method=TRUE)

model4.4.14_Laf<-glmmTMB(log_FH_Laf~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Laf_FH_NoZero, na.action=na.omit)
summary(model4.4.14_Laf)
Anova(model4.4.14_Laf, p.adjust.method=TRUE)

model4.5.14_Laf<-glmmTMB(sqrt_SS_Laf~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Laf_SS_NoZero, na.action=na.omit)
summary(model4.5.14_Laf)
Anova(model4.5.14_Laf, p.adjust.method=TRUE)

model4.6.14_Laf<-glmmTMB(log_MA_Laf~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.6.14_Laf)
Anova(model4.6.14_Laf, p.adjust.method=TRUE)

model4.7.14_Laf<-glmmTMB(log_GR_Laf~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Laf, na.action=na.omit)
summary(model4.7.14_Laf)
Anova(model4.7.14_Laf, p.adjust.method=TRUE)

#***************************** 
# Q4c. Uppsala
#*****************************

# Subset full experiment data to just data from the Uppsala common garden
CG_Upp=subset(CG_Data, Garden=="Uppsala")

# Subset CG_Upp to analyze zero-inflated flower head numbers as above
CG_Upp_FH_NoZero=subset(CG_Upp, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# Subset CG_Upp to analyze zero-inflated seed set masses as above
CG_Upp_SS_NoZero=subset(CG_Upp, !SeedSetMassFull=="0") |> 
  droplevels()

# Perform appropriate transformations for all quantitative variables below
# Each quantitative variable is also used in its z-score standardized form for equal comparisons across all models
sqrt_FH_Upp<-sqrt(CG_Upp_FH_NoZero$FH_Full_z)
log_SS_Upp<-log(CG_Upp_SS_NoZero$SS_Full_z)
sqrt_MA_Upp<-sqrt(CG_Upp$MA_1st_z)

# Q4c. BIO1
model4.1.1_Upp<-glmmTMB(Survived1st~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.1.1_Upp)
Anova(model4.1.1_Upp, p.adjust.method=TRUE)

model4.2.1_Upp<-glmmTMB(FloweredFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.2.1_Upp)
Anova(model4.2.1_Upp, p.adjust.method=TRUE)

model4.3.1_Upp<-glmmTMB(SeededFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.3.1_Upp)
Anova(model4.3.1_Upp, p.adjust.method=TRUE)

model4.4.1_Upp<-glmmTMB(sqrt_FH_Upp~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Upp_FH_NoZero, na.action=na.omit)
summary(model4.4.1_Upp)
Anova(model4.4.1_Upp, p.adjust.method=TRUE)

model4.5.1_Upp<-glmmTMB(log_SS_Upp~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Upp_SS_NoZero, na.action=na.omit)
summary(model4.5.1_Upp)
Anova(model4.5.1_Upp, p.adjust.method=TRUE)

model4.6.1_Upp<-glmmTMB(sqrt_MA_Upp~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.6.1_Upp)
Anova(model4.6.1_Upp, p.adjust.method=TRUE)

model4.7.1_Upp<-glmmTMB(GR_1st_z~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.7.1_Upp)
Anova(model4.7.1_Upp, p.adjust.method=TRUE)

# Q4c. BIO5
model4.1.5_Upp<-glmmTMB(Survived1st~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.1.5_Upp)
Anova(model4.1.5_Upp, p.adjust.method=TRUE)

model4.2.5_Upp<-glmmTMB(FloweredFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.2.5_Upp)
Anova(model4.2.5_Upp, p.adjust.method=TRUE)

model4.3.5_Upp<-glmmTMB(SeededFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.3.5_Upp)
Anova(model4.3.5_Upp, p.adjust.method=TRUE)

model4.4.5_Upp<-glmmTMB(sqrt_FH_Upp~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Upp_FH_NoZero, na.action=na.omit)
summary(model4.4.5_Upp)
Anova(model4.4.5_Upp, p.adjust.method=TRUE)

model4.5.5_Upp<-glmmTMB(log_SS_Upp~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Upp_SS_NoZero, na.action=na.omit)
summary(model4.5.5_Upp)
Anova(model4.5.5_Upp, p.adjust.method=TRUE)

model4.6.5_Upp<-glmmTMB(sqrt_MA_Upp~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.6.5_Upp)
Anova(model4.6.5_Upp, p.adjust.method=TRUE)

model4.7.5_Upp<-glmmTMB(GR_1st_z~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.7.5_Upp)
Anova(model4.7.5_Upp, p.adjust.method=TRUE)

# Q4c. BIO6
model4.1.6_Upp<-glmmTMB(Survived1st~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.1.6_Upp)
Anova(model4.1.6_Upp, p.adjust.method=TRUE)

model4.2.6_Upp<-glmmTMB(FloweredFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.2.6_Upp)
Anova(model4.2.6_Upp, p.adjust.method=TRUE)

model4.3.6_Upp<-glmmTMB(SeededFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.3.6_Upp)
Anova(model4.3.6_Upp, p.adjust.method=TRUE)

model4.4.6_Upp<-glmmTMB(sqrt_FH_Upp~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Upp_FH_NoZero, na.action=na.omit)
summary(model4.4.6_Upp)
Anova(model4.4.6_Upp, p.adjust.method=TRUE)

model4.5.6_Upp<-glmmTMB(log_SS_Upp~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Upp_SS_NoZero, na.action=na.omit)
summary(model4.5.6_Upp)
Anova(model4.5.6_Upp, p.adjust.method=TRUE)

model4.6.6_Upp<-glmmTMB(sqrt_MA_Upp~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.6.6_Upp)
Anova(model4.6.6_Upp, p.adjust.method=TRUE)

model4.7.6_Upp<-glmmTMB(GR_1st_z~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.7.6_Upp)
Anova(model4.7.6_Upp, p.adjust.method=TRUE)

# Q4c. BIO12
model4.1.12_Upp<-glmmTMB(Survived1st~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.1.12_Upp)
Anova(model4.1.12_Upp, p.adjust.method=TRUE)

model4.2.12_Upp<-glmmTMB(FloweredFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.2.12_Upp)
Anova(model4.2.12_Upp, p.adjust.method=TRUE)

model4.3.12_Upp<-glmmTMB(SeededFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.3.12_Upp)
Anova(model4.3.12_Upp, p.adjust.method=TRUE)

model4.4.12_Upp<-glmmTMB(sqrt_FH_Upp~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Upp_FH_NoZero, na.action=na.omit)
summary(model4.4.12_Upp)
Anova(model4.4.12_Upp, p.adjust.method=TRUE)

model4.5.12_Upp<-glmmTMB(log_SS_Upp~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Upp_SS_NoZero, na.action=na.omit)
summary(model4.5.12_Upp)
Anova(model4.5.12_Upp, p.adjust.method=TRUE)

model4.6.12_Upp<-glmmTMB(sqrt_MA_Upp~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.6.12_Upp)
Anova(model4.6.12_Upp, p.adjust.method=TRUE)

model4.7.12_Upp<-glmmTMB(GR_1st_z~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.7.12_Upp)
Anova(model4.7.12_Upp, p.adjust.method=TRUE)

# Q4c. BIO13
model4.1.13_Upp<-glmmTMB(Survived1st~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.1.13_Upp)
Anova(model4.1.13_Upp, p.adjust.method=TRUE)

model4.2.13_Upp<-glmmTMB(FloweredFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.2.13_Upp)
Anova(model4.2.13_Upp, p.adjust.method=TRUE)

model4.3.13_Upp<-glmmTMB(SeededFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.3.13_Upp)
Anova(model4.3.13_Upp, p.adjust.method=TRUE)

model4.4.13_Upp<-glmmTMB(sqrt_FH_Upp~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Upp_FH_NoZero, na.action=na.omit)
summary(model4.4.13_Upp)
Anova(model4.4.13_Upp, p.adjust.method=TRUE)

model4.5.13_Upp<-glmmTMB(log_SS_Upp~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Upp_SS_NoZero, na.action=na.omit)
summary(model4.5.13_Upp)
Anova(model4.5.13_Upp, p.adjust.method=TRUE)

model4.6.13_Upp<-glmmTMB(sqrt_MA_Upp~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.6.13_Upp)
Anova(model4.6.13_Upp, p.adjust.method=TRUE)

model4.7.13_Upp<-glmmTMB(GR_1st_z~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.7.13_Upp)
Anova(model4.7.13_Upp, p.adjust.method=TRUE)

# Q4c. BIO14
model4.1.14_Upp<-glmmTMB(Survived1st~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.1.14_Upp)
Anova(model4.1.14_Upp, p.adjust.method=TRUE)

model4.2.14_Upp<-glmmTMB(FloweredFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.2.14_Upp)
Anova(model4.2.14_Upp, p.adjust.method=TRUE)

model4.3.14_Upp<-glmmTMB(SeededFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Upp, family = "binomial", na.action=na.omit)
summary(model4.3.14_Upp)
Anova(model4.3.14_Upp, p.adjust.method=TRUE)

model4.4.14_Upp<-glmmTMB(sqrt_FH_Upp~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Upp_FH_NoZero, na.action=na.omit)
summary(model4.4.14_Upp)
Anova(model4.4.14_Upp, p.adjust.method=TRUE)

model4.5.14_Upp<-glmmTMB(log_SS_Upp~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Upp_SS_NoZero, na.action=na.omit)
summary(model4.5.14_Upp)
Anova(model4.5.14_Upp, p.adjust.method=TRUE)

model4.6.14_Upp<-glmmTMB(sqrt_MA_Upp~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.6.14_Upp)
Anova(model4.6.14_Upp, p.adjust.method=TRUE)

model4.7.14_Upp<-glmmTMB(GR_1st_z~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Upp, na.action=na.omit)
summary(model4.7.14_Upp)
Anova(model4.7.14_Upp, p.adjust.method=TRUE)

#***************************** 
# Q4d. Montpellier
#*****************************

# Subset full experiment data to just data from the Montpellier common garden
CG_Mon=subset(CG_Data, Garden=="Montpellier")

# Subset CG_Mon to analyze zero-inflated flower head numbers as above
CG_Mon_FH_NoZero=subset(CG_Mon, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# Subset CG_Mon to analyze zero-inflated seed set masses as above
CG_Mon_SS_NoZero=subset(CG_Mon, !SeedSetMassFull=="0") |> 
  droplevels()

# Perform appropriate transformations for all quantitative variables below
# Each quantitative variable is also used in its z-score standardized form for equal comparisons across all models
sqrt_FH_Mon<-sqrt(CG_Mon_FH_NoZero$FH_Full_z)
log_SS_Mon<-log(CG_Mon_SS_NoZero$SS_Full_z)
sqrt_MA_Mon<-sqrt(CG_Mon$MA_1st_z)

# Q4d. BIO1
model4.1.1_Mon<-glmmTMB(Survived1st~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.1.1_Mon)
Anova(model4.1.1_Mon, p.adjust.method=TRUE)

model4.2.1_Mon<-glmmTMB(FloweredFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.2.1_Mon)
Anova(model4.2.1_Mon, p.adjust.method=TRUE)

model4.3.1_Mon<-glmmTMB(SeededFull~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.3.1_Mon)
Anova(model4.3.1_Mon, p.adjust.method=TRUE)

model4.4.1_Mon<-glmmTMB(sqrt_FH_Mon~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mon_FH_NoZero, na.action=na.omit)
summary(model4.4.1_Mon)
Anova(model4.4.1_Mon, p.adjust.method=TRUE)

model4.5.1_Mon<-glmmTMB(log_SS_Mon~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mon_SS_NoZero, na.action=na.omit)
summary(model4.5.1_Mon)
Anova(model4.5.1_Mon, p.adjust.method=TRUE)

model4.6.1_Mon<-glmmTMB(sqrt_MA_Mon~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.6.1_Mon)
Anova(model4.6.1_Mon, p.adjust.method=TRUE)

model4.7.1_Mon<-glmmTMB(GR_1st_z~bio1diff+I(bio1diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.7.1_Mon)
Anova(model4.7.1_Mon, p.adjust.method=TRUE)

# Q4d. BIO5
model4.1.5_Mon<-glmmTMB(Survived1st~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.1.5_Mon)
Anova(model4.1.5_Mon, p.adjust.method=TRUE)

model4.2.5_Mon<-glmmTMB(FloweredFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.2.5_Mon)
Anova(model4.2.5_Mon, p.adjust.method=TRUE)

model4.3.5_Mon<-glmmTMB(SeededFull~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.3.5_Mon)
Anova(model4.3.5_Mon, p.adjust.method=TRUE)

model4.4.5_Mon<-glmmTMB(sqrt_FH_Mon~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mon_FH_NoZero, na.action=na.omit)
summary(model4.4.5_Mon)
Anova(model4.4.5_Mon, p.adjust.method=TRUE)

model4.5.5_Mon<-glmmTMB(log_SS_Mon~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mon_SS_NoZero, na.action=na.omit)
summary(model4.5.5_Mon)
Anova(model4.5.5_Mon, p.adjust.method=TRUE)

model4.6.5_Mon<-glmmTMB(sqrt_MA_Mon~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.6.5_Mon)
Anova(model4.6.5_Mon, p.adjust.method=TRUE)

model4.7.5_Mon<-glmmTMB(GR_1st_z~bio5diff+I(bio5diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.7.5_Mon)
Anova(model4.7.5_Mon, p.adjust.method=TRUE)

# Q4d. BIO6
model4.1.6_Mon<-glmmTMB(Survived1st~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.1.6_Mon)
Anova(model4.1.6_Mon, p.adjust.method=TRUE)

model4.2.6_Mon<-glmmTMB(FloweredFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.2.6_Mon)
Anova(model4.2.6_Mon, p.adjust.method=TRUE)

model4.3.6_Mon<-glmmTMB(SeededFull~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.3.6_Mon)
Anova(model4.3.6_Mon, p.adjust.method=TRUE)

model4.4.6_Mon<-glmmTMB(sqrt_FH_Mon~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mon_FH_NoZero, na.action=na.omit)
summary(model4.4.6_Mon)
Anova(model4.4.6_Mon, p.adjust.method=TRUE)

model4.5.6_Mon<-glmmTMB(log_SS_Mon~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mon_SS_NoZero, na.action=na.omit)
summary(model4.5.6_Mon)
Anova(model4.5.6_Mon, p.adjust.method=TRUE)

model4.6.6_Mon<-glmmTMB(sqrt_MA_Mon~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.6.6_Mon)
Anova(model4.6.6_Mon, p.adjust.method=TRUE)

model4.7.6_Mon<-glmmTMB(GR_1st_z~bio6diff+I(bio6diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.7.6_Mon)
Anova(model4.7.6_Mon, p.adjust.method=TRUE)

# Q4d. BIO12
model4.1.12_Mon<-glmmTMB(Survived1st~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.1.12_Mon)
Anova(model4.1.12_Mon, p.adjust.method=TRUE)

model4.2.12_Mon<-glmmTMB(FloweredFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.2.12_Mon)
Anova(model4.2.12_Mon, p.adjust.method=TRUE)

model4.3.12_Mon<-glmmTMB(SeededFull~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.3.12_Mon)
Anova(model4.3.12_Mon, p.adjust.method=TRUE)

model4.4.12_Mon<-glmmTMB(sqrt_FH_Mon~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mon_FH_NoZero, na.action=na.omit)
summary(model4.4.12_Mon)
Anova(model4.4.12_Mon, p.adjust.method=TRUE)

model4.5.12_Mon<-glmmTMB(log_SS_Mon~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mon_SS_NoZero, na.action=na.omit)
summary(model4.5.12_Mon)
Anova(model4.5.12_Mon, p.adjust.method=TRUE)

model4.6.12_Mon<-glmmTMB(sqrt_MA_Mon~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.6.12_Mon)
Anova(model4.6.12_Mon, p.adjust.method=TRUE)

model4.7.12_Mon<-glmmTMB(GR_1st_z~bio12diff+I(bio12diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.7.12_Mon)
Anova(model4.7.12_Mon, p.adjust.method=TRUE)

# Q4d. BIO13
model4.1.13_Mon<-glmmTMB(Survived1st~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.1.13_Mon)
Anova(model4.1.13_Mon, p.adjust.method=TRUE)

model4.2.13_Mon<-glmmTMB(FloweredFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.2.13_Mon)
Anova(model4.2.13_Mon, p.adjust.method=TRUE)

model4.3.13_Mon<-glmmTMB(SeededFull~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.3.13_Mon)
Anova(model4.3.13_Mon, p.adjust.method=TRUE)

model4.4.13_Mon<-glmmTMB(sqrt_FH_Mon~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mon_FH_NoZero, na.action=na.omit)
summary(model4.4.13_Mon)
Anova(model4.4.13_Mon, p.adjust.method=TRUE)

model4.5.13_Mon<-glmmTMB(log_SS_Mon~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mon_SS_NoZero, na.action=na.omit)
summary(model4.5.13_Mon)
Anova(model4.5.13_Mon, p.adjust.method=TRUE)

model4.6.13_Mon<-glmmTMB(sqrt_MA_Mon~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.6.13_Mon)
Anova(model4.6.13_Mon, p.adjust.method=TRUE)

model4.7.13_Mon<-glmmTMB(GR_1st_z~bio13diff+I(bio13diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.7.13_Mon)
Anova(model4.7.13_Mon, p.adjust.method=TRUE)

# Q4d. BIO14
model4.1.14_Mon<-glmmTMB(Survived1st~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.1.14_Mon)
Anova(model4.1.14_Mon, p.adjust.method=TRUE)

model4.2.14_Mon<-glmmTMB(FloweredFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.2.14_Mon)
Anova(model4.2.14_Mon, p.adjust.method=TRUE)

model4.3.14_Mon<-glmmTMB(SeededFull~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mon, family = "binomial", na.action=na.omit)
summary(model4.3.14_Mon)
Anova(model4.3.14_Mon, p.adjust.method=TRUE)

model4.4.14_Mon<-glmmTMB(sqrt_FH_Mon~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mon_FH_NoZero, na.action=na.omit)
summary(model4.4.14_Mon)
Anova(model4.4.14_Mon, p.adjust.method=TRUE)

model4.5.14_Mon<-glmmTMB(log_SS_Mon~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mon_SS_NoZero, na.action=na.omit)
summary(model4.5.14_Mon)
Anova(model4.5.14_Mon, p.adjust.method=TRUE)

model4.6.14_Mon<-glmmTMB(sqrt_MA_Mon~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.6.14_Mon)
Anova(model4.6.14_Mon, p.adjust.method=TRUE)

model4.7.14_Mon<-glmmTMB(GR_1st_z~bio14diff+I(bio14diff^2)+(1|Population), data=CG_Mon, na.action=na.omit)
summary(model4.7.14_Mon)
Anova(model4.7.14_Mon, p.adjust.method=TRUE)


#***************************** 
# 2. EXPERIMENT DATA ANALYSIS - e) MODELS TO ADDRESS Q5
#*****************************

# Now analyze limited dataset for spatial lags in adaptation to climate change

# Step 1: # Data will be first analyzed using both linear and quadratic terms for bio1 (Mean Annual Temperature, raw values)
    # Run glmmTMBs for each of the 7 response variables in the order listed above (excluding Herbivory)
    # Each of the 6 bioclimatic variables (both linear and quadratic terms incuded as predictors) run in a separate model with the number system listed:
      # 1 = Mean Annual Temperature
      # 5 = Max Temperature of Warmest Month
      # 6 = Min Temperature of Coldest Month
      # 12 = Annual Precipitation
      # 13 = Precipitation of Wettest Month
      # 14 = Precipitation of Driest Month
    # Population included as a random effect
# Step 2: Run model summaries and Anova to obtain test statistics and p values
# Step 3: Models will then be repeated without the linear term (NL) to test for a lateral shift in the quadratic curve using anova and AIC
# Step 4: Plot data for each model with a significant linear or quadratic relationship

# Subset data to only the North American populations in the Mississauga garden
CG_Mis_NA=subset(CG_Mis, ContOrigin=="NAM") |> 
  droplevels()

# Subset CG_Mis_NA to analyze zero-inflated flower head numbers as above
CG_Mis_NA_FH_NoZero=subset(CG_Mis_NA, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# Subset CG_Mis_NA analyze zero-inflated seed set masses as above
CG_Mis_NA_SS_NoZero=subset(CG_Mis_NA, !SeedSetMassFull=="0") |> 
  droplevels()

# Perform log and sqrt transformations for all quantitative variables below (not all will be used but all tested)
log_FH_Mis_NA<-log(CG_Mis_NA_FH_NoZero$FlowerHeadNumberFull)
log_SS_Mis_NA<-log(CG_Mis_NA_SS_NoZero$SeedSetMassFull)
log_MA_Mis_NA<-log(CG_Mis_NA$MaxArea1st)

# Full model for Survived1st
model5.1_Mis<-glmmTMB(Survived1st~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Mis_NA, family = "binomial", na.action=na.omit)
summary(model5.1_Mis)
Anova(model5.1_Mis, p.adjust.method=TRUE)

# Reduced model (no linear term) for Survived 1st
model5.1_MisNL<-glmmTMB(Survived1st~I(MATdiff^2)+(1|Population), data=CG_Mis_NA, family = "binomial", na.action=na.omit)
summary(model5.1_MisNL)

anova(model5.1_Mis, model5.1_MisNL) # anova indicates significant difference between the models
AIC(model5.1_Mis)
AIC(model5.1_MisNL) # > AIC(model5.1_Mis) indicates full model is the better fit than reduced model
# anova and AIC values indicate significant shift away from the origin of the quadratic curve

Fig5.1<-CG_Mis_NA |>
  group_by(MATdiff,Population) |>
  dplyr::summarize(Surv.mean = mean(Survived1st)) |>
  ggplot(aes(x = MATdiff,  y = Surv.mean, color = "black")) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", formula = y ~ x + I(x^2), se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    formula = y ~ x + I(x^2),
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax..)) +
  ggtitle("A") +
  ylab("Proportion Survived") +
  xlab("Δ Mean Annual Temperature (°C)") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette1, name = "") + 
  scale_fill_manual(values = Palette1, name = "") + 
  coord_fixed(ratio = 19.230769231) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig5.1 <- Fig5.1 + theme(legend.position = "none")
Fig5.1

# Full model for FloweredFull
model5.2_Mis<-glmmTMB(FloweredFull~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Mis_NA, family = "binomial", na.action=na.omit)
summary(model5.2_Mis)
Anova(model5.2_Mis, p.adjust.method=TRUE)

Fig5.2<-CG_Mis_NA |>
  group_by(MATdiff, Population) |>
  dplyr::summarize(Flowered.mean = mean(FloweredFull)) |>
  ggplot(aes(x = MATdiff,  y = Flowered.mean, color = "black")) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax..)) +
  ggtitle("B") +
  ylab("Proportion Flowered") +
  xlab("Δ Mean Annual Temperature (°C)") +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette1, name = "") + 
  scale_fill_manual(values = Palette1, name = "") + 
  coord_fixed(ratio = 19.230769231) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig5.2 <- Fig5.2 + theme(legend.position = "none")
Fig5.2

# Full model for SeededFull
model5.3_Mis<-glmmTMB(SeededFull~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Mis_NA, family = "binomial", na.action=na.omit)
summary(model5.3_Mis)
Anova(model5.3_Mis, p.adjust.method=TRUE)

Fig5.3<-CG_Mis_NA |>
  group_by(MATdiff, Population) |>
  dplyr::summarize(Seeded.mean = mean(SeededFull)) |>
  ggplot(aes(x = MATdiff,  y = Seeded.mean, color = "black")) +
  geom_point(size = 1.5) +
  stat_smooth(method = "glm", se = FALSE, method.args=list(family="binomial"), fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "glm",
    se = TRUE,
    method.args=list(family="binomial"),
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.02), ymax = ..ymax..)) +
  ggtitle("C") +
  ylab("Proportion Producing Seeds") +
  xlab("Δ Mean Annual Temperature (°C)") +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_y_continuous(limits = c(-0.02,1.02), breaks = c(0,0.2,0.4,0.6,0.8,1.0), expand = c(0,0)) +
  scale_color_manual(values = Palette1, name = "") + 
  scale_fill_manual(values = Palette1, name = "") + 
  coord_fixed(ratio = 19.230769231) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig5.3 <- Fig5.3 + theme(legend.position = "none")
Fig5.3

# Full model for log_FH_Mis_NA
model5.4_Mis<-glmmTMB(log_FH_Mis_NA~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Mis_NA_FH_NoZero, na.action=na.omit)
summary(model5.4_Mis)
Anova(model5.4_Mis, p.adjust.method=TRUE)

# Full model for log_SS_Mis_NA
model5.5_Mis<-glmmTMB(log_SS_Mis_NA~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Mis_NA_SS_NoZero, na.action=na.omit)
summary(model5.5_Mis)
Anova(model5.5_Mis, p.adjust.method=TRUE)

# Reduced model (no linear term) for log_SS_Mis_NA
model5.5_MisNL<-glmmTMB(log_SS_Mis_NA~I(MATdiff^2)+(1|Population), data=CG_Mis_NA_SS_NoZero, na.action=na.omit)
summary(model5.5_MisNL)

anova(model5.5_Mis, model5.5_MisNL) # anova indicates non-significant difference between the models
AIC(model5.5_Mis)
AIC(model5.5_MisNL) # < AIC(model5.5_Mis) indicates reduced model is the better fit than full model
# anova and AIC values indicate non-significant shift away from the origin of the quadratic curve

Fig5.5<-CG_Mis_NA_SS_NoZero |>
  filter(!is.na(SeedSetMassFull)) |>
  group_by(MATdiff, Population) |>
  dplyr::summarize(SS.mean = mean(SeedSetMassFull)) |>
  ggplot(aes(x = MATdiff,  y = SS.mean, color = "black")) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -0.2), ymax = ..ymax..)) +
  ggtitle("D") +
  ylab("Seed Set Mass (g)") +
  xlab("Δ Mean Annual Temperature (°C)") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_y_continuous(limits = c(-0.2,10.2), breaks = c(0,2,4,6,8,10), expand = c(0,0)) +
  scale_color_manual(values = Palette1, name = "") + 
  scale_fill_manual(values = Palette1, name = "") + 
  coord_fixed(ratio = 1.9230769231) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig5.5 <- Fig5.5 + theme(legend.position = "none")
Fig5.5

# Full model for log_MA_Mis_NA
model5.6_Mis<-glmmTMB(log_MA_Mis_NA~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Mis_NA, na.action=na.omit)
summary(model5.6_Mis)
Anova(model5.6_Mis, p.adjust.method=TRUE)

# Reduced model (no linear term) for log_MA_Mis_NA
model5.6_MisNL<-glmmTMB(log_MA_Mis_NA~I(MATdiff^2)+(1|Population), data=CG_Mis_NA, na.action=na.omit)
summary(model5.6_MisNL)

anova(model5.6_Mis, model5.6_MisNL)  # anova indicates marginally significant difference between the models
AIC(model5.6_Mis)
AIC(model5.6_MisNL) # < AIC(model5.5_Mis) indicates full model is the better fit than reduced model
# anova and AIC values indicate non-significant shift away from the origin of the quadratic curve

Fig5.6<-CG_Mis_NA |>
  filter(!is.na(MaxArea1st)) |>
  group_by(MATdiff, Population) |>
  dplyr::summarize(MA.mean = mean(MaxArea1st)) |>
  ggplot(aes(x = MATdiff,  y = MA.mean, color = "black")) +
  geom_point(size = 1.5) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  stat_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE,
    fullrange = TRUE,
    geom = "ribbon",
    alpha = 0.3,
    color = NA,
    aes(ymin = pmax(..ymin.., -6), ymax = ..ymax..)) +
  ggtitle("E") +
  ylab("Maximum Plant Area (cm²)") +
  xlab("Δ Mean Annual Temperature (°C)") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_y_continuous(limits = c(-6,306), breaks = c(0,100,200,300), expand = c(0,0)) +
  scale_color_manual(values = Palette1, name = "") + 
  scale_fill_manual(values = Palette1, name = "") + 
  coord_fixed(ratio = 0.0641025641) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=12))
Fig5.6 <- Fig5.6 + theme(legend.position = "none")
Fig5.6

# Full model for GrowthRate1st
model5.7_Mis<-glmmTMB(GrowthRate1st~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Mis_NA, na.action=na.omit)
summary(model5.7_Mis)
Anova(model5.7_Mis, p.adjust.method=TRUE)

# Now repeat for Uppsala
# Subset data to only the European populations in the Uppsala garden
CG_Upp_EU=subset(CG_Upp, ContOrigin=="EUR") |> 
  droplevels()

# Subset CG_Upp_EU to analyze zero-inflated flower head numbers as above
CG_Upp_EU_FH_NoZero=subset(CG_Upp_EU, !FlowerHeadNumberFull=="0") |> 
  droplevels()

# Subset CG_Upp_EU analyze zero-inflated seed set masses as above
CG_Upp_EU_SS_NoZero=subset(CG_Upp_EU, !SeedSetMassFull=="0") |> 
  droplevels()

# Perform log and sqrt transformations for all quantitative variables below (not all will be used but all tested)
log_FH_Upp_EU<-log(CG_Upp_EU_FH_NoZero$FlowerHeadNumberFull)
log_SS_Upp_EU<-log(CG_Upp_EU_SS_NoZero$SeedSetMassFull)
sqrt_MA_Upp_EU<-sqrt(CG_Upp_EU$MaxArea1st)

model5.1_Upp<-glmmTMB(Survived1st~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Upp_EU, family = "binomial", na.action=na.omit)
summary(model5.1_Upp)
Anova(model5.1_Upp, p.adjust.method=TRUE)

model5.2_Upp<-glmmTMB(FloweredFull~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Upp_EU, family = "binomial", na.action=na.omit)
summary(model5.2_Upp)
Anova(model5.2_Upp, p.adjust.method=TRUE)

model5.3_Upp<-glmmTMB(SeededFull~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Upp_EU, family = "binomial", na.action=na.omit)
summary(model5.3_Upp)
Anova(model5.3_Upp, p.adjust.method=TRUE)

model5.4_Upp<-glmmTMB(log_FH_Upp_EU~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Upp_EU_FH_NoZero, na.action=na.omit)
summary(model5.4_Upp)
Anova(model5.4_Upp, p.adjust.method=TRUE)

model5.5_Upp<-glmmTMB(log_SS_Upp_EU~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Upp_EU_SS_NoZero, na.action=na.omit)
summary(model5.5_Upp)
Anova(model5.5_Upp, p.adjust.method=TRUE)

model5.6_Upp<-glmmTMB(sqrt_MA_Upp_EU~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Upp_EU, na.action=na.omit)
summary(model5.6_Upp)
Anova(model5.6_Upp, p.adjust.method=TRUE)

model5.7_Upp<-glmmTMB(GrowthRate1st~MATdiff+I(MATdiff^2)+(1|Population), data=CG_Upp_EU, na.action=na.omit)
summary(model5.7_Upp)
Anova(model5.7_Upp, p.adjust.method=TRUE)

#***************************** 
# 3. SUPPLEMENTAL ANALYSIS - Climatic distance (CD) model
#*****************************

# Run linear model to test average CDs among sampled populations within each continent transplanted into common gardens on each continent
CDmodel <- lm(WeightedCD ~ ContOrigin*ContGarden, data=CG_Data, na.action=na.omit)
summary(CDmodel)
Anova(CDmodel, p.adjust.method=TRUE)

CDmodel.em<-emmeans(CDmodel, ~ ContOrigin*ContGarden)
print(CDmodel.em)
multcomp::cld(CDmodel.em)
# European populations in Europe experience the smallest average CD between environment of origin and common garden environment
