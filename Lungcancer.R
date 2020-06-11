#Loading libraries 
library(RTCGA)
library(RTCGA.clinical)
library(survival)
library(survminer)
library(gridExtra)
library(kableExtra)
library(tigerstats)
library(tidyverse)
library(broom)

#extracting data from TCGA 
LUSC <- survivalTCGA(LUSC.clinical,extract.cols = c("patient.gender","patient.ethnicity","patient.race","patient.age_at_initial_pathologic_diagnosis",
                                                    "patient.stage_event.pathologic_stage", "patient.stage_event.tnm_categories.pathologic_categories.pathologic_m",
                                                    "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n","patient.stage_event.tnm_categories.pathologic_categories.pathologic_t",
                                                    "patient.tobacco_smoking_history"))

LUAD <- survivalTCGA(LUAD.clinical,extract.cols= c("patient.gender","patient.ethnicity","patient.race","patient.age_at_initial_pathologic_diagnosis",
                                                   "patient.stage_event.pathologic_stage", "patient.stage_event.tnm_categories.pathologic_categories.pathologic_m",
                                                   "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n","patient.stage_event.tnm_categories.pathologic_categories.pathologic_t",
                                                   "patient.tobacco_smoking_history"))

# A thir dataframe which consists of 4 columns: times,patient_code,vital_status, and disease code(LUSC/LUAD)
LUSC_LUAD <-survivalTCGA(LUSC.clinical,LUAD.clinical,extract.cols= c("admin.disease_code"))


#Data Wrangling
## Step1 - Eliminate negative values from the times column and convert values from days to years

LUSC_modified <- LUSC %>% filter(times > 0) %>% mutate(times=times/365)
LUAD_modified <- LUAD %>% filter(times > 0) %>% mutate(times=times/365)
LUSC_LUAD_modified <- LUSC_LUAD %>% filter(times > 0) %>% mutate(times=times/365)
 

#modify column names
#The column names are too long. So, I will modify them for ease of usage 

colnames(LUSC_modified) <- c("time_to_event","ID","vital_status","gender","ethnicity","race","age","pathologic_stage","stage_m","stage_n","stage_t","smoking_history")
colnames(LUAD_modified) <- c("time_to_event","ID","vital_status","gender","ethnicity","race","age","pathologic_stage","stage_m","stage_n","stage_t","smoking_history")


## Step2 - Check for missing data
LUSC_missing <- LUSC_modified %>% sapply(function(x) round((sum(is.na(x))/nrow(LUSC_modified))*100,digits=2))
LUSC_missing <- as.data.frame(LUSC_missing)
LUSC_missing <- cbind(variable=rownames(LUSC_missing),data.frame(LUSC_missing,row.names=NULL))

LUAD_missing <- LUAD_modified %>% sapply(function(x) round((sum(is.na(x))/nrow(LUAD_modified))*100,digits=2))
LUAD_missing <- as.data.frame(LUAD_missing)

#combining the dataframes to visualize the data
LUSC_LUAD_missing <- cbind(LUSC_missing,LUAD_missing,row.names=NULL)
gather(LUSC_LUAD_missing,type,prop,LUSC_missing:LUAD_missing) %>% filter(prop>0.00) %>% ggplot( aes(x=variable,y=prop,fill=type))+geom_bar(stat="identity",position="dodge")+coord_flip()+theme_classic()+xlab("")+ylab("% missing")+
scale_fill_discrete(name="Cancer",labels=c("LUSC","LUAD"))+ggtitle("Missing values in LUSC and LUAD")


# For race and ethinicity variables, LUSC has 10% and 20% missing values while LUAD has 20% and 34% missing values respectively 
# Let's explore the RACE and ETHNICITY variables** 
ethnicity_LUSC_plot <- LUSC_modified %>% arrange(ethnicity) %>% mutate(index=1:n()) %>% 
ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=ethnicity,shape=factor(vital_status)))+geom_segment()+geom_point()+
xlab("Time to Event(years)")+ylab("")+scale_shape_discrete(name="Status",labels=c("Censored","Death"))+theme_classic()

ethnicity_LUAD_plot <- LUAD_modified %>% arrange(ethnicity) %>% mutate(index=1:n()) %>% 
ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=ethnicity,shape=factor(vital_status)))+geom_segment()+geom_point()+
xlab("Time to Event(years)")+ylab("")+theme(legend.position="none")+theme_classic()

#Explore race column 
race_LUSC_plot <- LUSC_modified %>% arrange(race) %>% mutate(index=1:n()) %>% 
ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=race,shape=factor(vital_status)))+geom_segment()+geom_point()+xlab("Time to Event(years)")+ylab("")+scale_color_discrete(name="Race(LUSC)")+scale_shape_discrete(name="Status",labels=c("Censored","Death"))+theme_classic()+labs(title="LUSC")

race_LUAD_plot <- LUAD_modified %>% arrange(race) %>% mutate(index=1:n()) %>% 
ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=race,shape=factor(vital_status)))+geom_segment()+geom_point()+
xlab("Time to Event(years)")+ylab("")+scale_color_discrete(name="Race(LUAD)")+scale_shape_discrete(name="Status",labels=c("Censored","Death"))+ theme_classic()+labs(title="LUAD")

race_LUSC_plot
kable(table(LUSC_modified$race),booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
race_LUAD_plot
kable(table(LUAD_modified$race),booktabs=T) %>% kable_styling(full_width=FALSE,position="left")

## Step3 - Data cleaning  

#Deleting rows with missing data in pathalogic_stage,stage_m,smoking_history
LUSC_nomissing_values <- LUSC_modified %>% select(!ethnicity) %>% filter(!is.na(pathologic_stage)) %>% filter(!is.na(stage_m)) %>% filter(!is.na(stage_n)) %>% filter(!is.na(smoking_history)) %>% filter(!is.na(race)) 
LUAD_nomissing_values <- LUAD_modified %>% select(!ethnicity) %>% filter(!is.na(pathologic_stage)) %>% filter(!is.na(stage_m)) %>% filter(!is.na(smoking_history)) %>% filter(!is.na(race)) 

#combining all the minority races together 
LUSC_nomissing_values <- LUSC_nomissing_values %>% mutate(race=ifelse(race!="white","other",race))
LUAD_nomissing_values <- LUAD_nomissing_values %>% mutate(race=ifelse(race!="white","other",race))

#Now that we have clean data, let's explore each variable 
#Data Exploration - for each category we will plot distributions
#Gender  
table_gender_LUSC <- xtabs(~gender+vital_status,data=LUSC_nomissing_values)
gender_LUSC <- LUSC_nomissing_values %>% arrange(gender) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=gender,shape=factor(vital_status)))+
  geom_segment()+geom_point()+ theme_classic()+xlab("Time to Event(years)")+ylab("")+  scale_shape_discrete(name="Status",labels=c("Censored","Event"))

table_gender_LUAD <- xtabs(~gender+vital_status,data=LUAD_nomissing_values)
gender_LUAD <- LUAD_nomissing_values %>% arrange(gender) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=gender,shape=factor(vital_status)))+
  geom_segment()+geom_point()+ theme_classic()+ xlab("Time to Event(years)")+ylab("")

gender_LUSC
gender_LUAD
#ggarrange(gender_LUSC,gender_LUAD,nrow=1,common.legend=TRUE,legend="bottom",labels=c("LUSC","LUAD"),hjust=0.0)

## Age
#Check age distribution
LUSC_nomissing_values$age <- as.numeric(LUSC_nomissing_values$age)
age_LUSC <- LUSC_nomissing_values %>% ggplot(aes(x=age))+geom_histogram(aes(y=..density..),color="black",fill="white")+
  geom_density(alpha=0.2,fill="#FF6666")+geom_vline(aes(xintercept=mean(age)),color="blue")+theme_classic()

LUAD_nomissing_values$age <- as.numeric(LUAD_nomissing_values$age)
age_LUAD<-LUAD_nomissing_values %>% ggplot(aes(x=age))+geom_histogram(aes(y=..density..),color="black",fill="white")+
  geom_density(alpha=0.2,fill="#FF6666")+geom_vline(aes(xintercept=mean(age)),color="blue")+theme_classic()

age_LUSC
age_LUAD
  
#Creating a new categorical variable agecat
LUSC_nomissing_values <- LUSC_nomissing_values %>% mutate(agecat=ifelse(age>=69,"above69","below69"))
LUAD_nomissing_values <- LUAD_nomissing_values %>% mutate(agecat=ifelse(age>=65,"above65","below65"))

#plotting categorical age distributions 
table_age_LUSC <- xtabs(~agecat+vital_status,data=LUSC_nomissing_values)
age_LUSC <- LUSC_nomissing_values %>% arrange(agecat) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=agecat,shape=factor(vital_status)))+
  geom_segment()+geom_point()+theme_classic()+xlab("Time to Event(years)")+ylab("")+  scale_shape_discrete(name="Status",labels=c("Censored","Event"))+theme(legend.title = element_blank())

table_age_LUAD <- xtabs(~agecat+vital_status,data=LUAD_nomissing_values)
age_LUAD <- LUAD_nomissing_values %>% arrange(agecat) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=agecat,shape=factor(vital_status)))+     geom_segment()+geom_point()+theme_classic()+xlab("Time to Event(years)")+ylab("")+guides(shape=FALSE)+theme(legend.title = element_blank())


#Pathologic stage
table_ps_LUSC<- ggtexttable(round(colPerc(xtabs(~pathologic_stage+vital_status,data=LUSC_nomissing_values)),digits = 0))
pathologic_stage_LUSC <- LUSC_nomissing_values %>% arrange(pathologic_stage) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=pathologic_stage,shape=factor(vital_status)))+
  geom_segment()+geom_point()+theme_classic()+xlab("Time to Event(years)")+ylab("")+
  scale_shape_discrete(name="Status",labels=c("Censored(0)","Event(1)"))+labs(title="LUSC")
pathologic_stage_LUSC <- ggarrange(pathologic_stage_LUSC,table_ps_LUSC,ncol=2,widths=c(2,1))
pathologic_stage_LUSC+annotate("text",x=0.85,y=0.85,label="% values for each stage")

table_ps_LUAD <- ggtexttable(round(colPerc(xtabs(~pathologic_stage+vital_status,data=LUAD_nomissing_values)),digits=0))
pathologic_stage_LUAD <- LUAD_nomissing_values %>% arrange(pathologic_stage) %>% mutate(index=1:n()) %>%     ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=pathologic_stage,shape=factor(vital_status)))+geom_segment()+geom_point()+theme_classic()+xlab("Time to Event(years)")+ylab("")+
  scale_shape_discrete(name="Status",labels=c("Censored (0)","Event(1)"))+labs(title="LUAD")

pathologic_stage_LUAD <- ggarrange(pathologic_stage_LUAD,table_ps_LUAD,ncol=2,widths=c(2,1))
pathologic_stage_LUAD+annotate("text",x=0.85,y=0.85,label="% values for each stage")

#Stage T,N,M
#LUSC 

#Stage_t
table_staget_LUSC <- tableGrob(as.table(round(colPerc(xtabs(~stage_t+vital_status,data=LUSC_nomissing_values)),digits=0)),theme=ttheme_default(base_size=10))
tstage_LUSC<- LUSC_nomissing_values %>% arrange(stage_t) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=stage_t,shape=factor(vital_status)))+
  geom_segment()+geom_point()+guides(shape=FALSE)+theme(legend.position = "top")+
  xlab("Time to Event(years)")+ylab("")+guides(shape=FALSE)+theme_classic()

#Stage_n
table_stagen_LUSC <- tableGrob(as.table(round(colPerc(xtabs(~stage_n+vital_status,data=LUSC_nomissing_values)),digits = 0)),theme=ttheme_default(base_size=10))
nstage_LUSC <- LUSC_nomissing_values %>% arrange(stage_n) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=stage_n,shape=factor(vital_status)))+
  geom_segment()+geom_point()+xlab("Time to Event(years)")+ylab("")+theme(legend.position = "top")+guides(shape=FALSE)+theme_classic()

# #Stage_m
table_stagem_LUSC<-tableGrob(as.table(round(colPerc(xtabs(~stage_m+vital_status,data=LUSC_nomissing_values)),digits=0)),theme=ttheme_default(base_size=10))
mstage_LUSC <- LUSC_nomissing_values %>% arrange(stage_m) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=stage_m,shape=factor(vital_status)))+ geom_segment()+geom_point()+
  xlab("Time to Event(years)")+ylab("")+theme(legend.position="bottom")+ scale_shape_discrete(name="Status",labels=c("Censored","Event"))+theme_classic()

tstage_LUSC <- ggarrange(tstage_LUSC,table_staget_LUSC)
tstage_LUSC+annotate("text",x=0.80,y=0.8,label="% values for each stage")
nstage_LUSC <- ggarrange(nstage_LUSC,table_stagen_LUSC)
nstage_LUSC+annotate("text",x=0.80,y=0.7,label="% values for each stage")
mstage_LUSC <- ggarrange(mstage_LUSC,table_stagem_LUSC)
mstage_LUSC+annotate("text",x=0.80,y=0.7,label="% values for each stage")


#LUAD
#Stage_t
table_staget_LUAD<- tableGrob(as.table(round(colPerc(xtabs(~stage_t+vital_status,data=LUAD_nomissing_values)),digits=0)),theme=ttheme_default(base_size=10))
tstage_LUAD<- LUAD_nomissing_values %>% arrange(stage_t) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=stage_t,shape=factor(vital_status)))+
  geom_segment()+geom_point()+guides(shape=FALSE)+theme(legend.position = "top")+
  xlab("Time to Event(years)")+ylab("")+guides(shape=FALSE)+theme_classic()

#Stage_n
table_stagen_LUAD <- tableGrob(as.table(round(colPerc(xtabs(~stage_n+vital_status,data=LUAD_nomissing_values)),digits = 0)),theme=ttheme_default(base_size=10))
nstage_LUAD <- LUAD_nomissing_values %>% arrange(stage_n) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=stage_n,shape=factor(vital_status)))+
  geom_segment()+geom_point()+xlab("Time to Event(years)")+ylab("")+theme(legend.position = "top")+guides(shape=FALSE)+theme_classic()

# #Stage_m
table_stagem_LUAD<-tableGrob(as.table(round(colPerc(xtabs(~stage_m+vital_status,data=LUAD_nomissing_values)),digits=0)),theme=ttheme_default(base_size=10))
mstage_LUAD <- LUAD_nomissing_values %>% arrange(stage_m) %>% mutate(index=1:n()) %>% ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=stage_m,shape=factor(vital_status)))+ geom_segment()+geom_point()+
  xlab("Time to Event(years)")+ylab("")+theme(legend.position="bottom")+ scale_shape_discrete(name="Status",labels=c("Censored","Event"))+theme_classic()

tstage_LUAD <- ggarrange(tstage_LUAD,table_staget_LUAD)
tstage_LUAD+annotate("text",x=0.80,y=0.8,label="% values for each stage")
nstage_LUAD <- ggarrange(nstage_LUAD,table_stagen_LUAD)
nstage_LUAD+annotate("text",x=0.80,y=0.7,label="% values for each stage")
mstage_LUAD <- ggarrange(mstage_LUAD,table_stagem_LUAD)
mstage_LUAD+annotate("text",x=0.80,y=0.7,label="% values for each stage")

#function to clean pathologic stages. Transforms subcategories to main category for cancer stage, stages t,n,m
change_pathologic_stages <- function(x){
  x <- x %>% mutate(pathologic_stage=ifelse(pathologic_stage=="stage ia","stage i",
                                            ifelse(pathologic_stage=="stage ib", "stage i",
                                                   ifelse(pathologic_stage =="stage iia", "stage ii",
                                                          ifelse(pathologic_stage=="stage iib", "stage ii",
                                                                 ifelse(pathologic_stage=="stage iiia", "stage iii",
                                                                        ifelse(pathologic_stage=="stage iiib", "stage iii", pathologic_stage))))))) 
  
  x<-  x %>% mutate(stage_m=ifelse(stage_m=="m1a","m1",
                                   ifelse(stage_m=="m1b","m1",stage_m)))
  
  x %>% mutate(stage_t=ifelse(stage_t=="t1a","t1",
                              ifelse(stage_t=="t1b", "t1",
                                     ifelse(stage_t=="t2a", "t2",
                                            ifelse(stage_t=="t2b", "t2",stage_t)))))
}

#Transforming pathologic stages
LUSC_nomissing_values <- change_pathologic_stages(LUSC_nomissing_values)
LUAD_nomissing_values <- change_pathologic_stages(LUAD_nomissing_values)

##Smoking History
table_sh_LUSC <- xtabs(~smoking_history+vital_status,data=LUSC_nomissing_values)
smoking_LUSC <- LUSC_nomissing_values %>% arrange(smoking_history) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=smoking_history,shape=factor(vital_status)))+ geom_segment()+geom_point()+scale_shape_discrete(name="Status",labels=c("Censored","Event"))+
  xlab("Time to Event(years)")+ylab("")+theme(legend.box = "vertical")+theme_classic()

table_sh_LUAD <-xtabs(~smoking_history+vital_status,data=LUAD_nomissing_values)
smoking_LUAD <- LUAD_nomissing_values %>% arrange(smoking_history) %>% mutate(index=1:n()) %>% 
  ggplot(aes(xend=0,y=index,x=time_to_event,yend=index,colour=smoking_history,shape=factor(vital_status)))+ geom_segment()+geom_point()+ xlab("Time to Event(years)")+ylab("")+theme_classic()+theme(legend.position = "none")

ggarrange(smoking_LUSC,smoking_LUAD,nrow=1,common.legend = TRUE,legend="bottom",labels=c("LUSC","LUAD"),hjust=0.0)
kable(table_sh_LUSC,"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
kable(table_sh_LUAD,"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="center")

#transforming the smoking history column
LUSC_nomissing_values <- LUSC_nomissing_values %>% mutate(smoking_history=ifelse(smoking_history=="current reformed smoker, duration not specified" | smoking_history=="current reformed smoker for > 15 years" | smoking_history=="current reformed smoker for < or = 15 years","reformed smoker",smoking_history))
LUAD_nomissing_values <- LUAD_nomissing_values %>% mutate(smoking_history=ifelse(smoking_history=="current reformed smoker, duration not specified" | smoking_history=="current reformed smoker for > 15 years" | smoking_history=="current reformed smoker for < or = 15 years","reformed smoker",smoking_history))

#converting the columns to factors
numeric <- c(1,6)
factor <- c(4,5,7,8,9,10,11,12)
LUSC_nomissing_values[,numeric] <- sapply(LUSC_nomissing_values[,numeric],as.numeric)
LUSC_nomissing_values[,factor] <- lapply(LUSC_nomissing_values[,factor],as.factor)
LUAD_nomissing_values[,numeric] <- sapply(LUAD_nomissing_values[,numeric],as.numeric)
LUAD_nomissing_values[,factor] <- lapply(LUAD_nomissing_values[,factor],as.factor)
LUSC_nomissing_values$smoking_history <- factor(LUSC_nomissing_values$smoking_history,levels=c("current smoker","reformed smoker","lifelong non-smoker"))
LUAD_nomissing_values$smoking_history <- factor(LUAD_nomissing_values$smoking_history,levels=c("current smoker","reformed smoker","lifelong non-smoker"))

#dropping values that do not have significant frequency in staget,n,m columns
LUSC_model <- LUSC_nomissing_values %>% filter(pathologic_stage!="stage iv") %>% droplevels()
LUSC_model <- LUSC_model %>% filter(stage_n!="n3") %>% droplevels()
LUSC_model <- LUSC_model %>% filter(stage_n!="nx") %>% droplevels()
LUSC_model <- LUSC_model %>% filter(stage_m!="m1") %>% droplevels()

LUAD_model <- LUAD_nomissing_values %>% filter(stage_t!="tx") %>% droplevels()
LUAD_model <- LUAD_model %>% filter(stage_n!="n3") %>% droplevels()
LUAD_model <- LUAD_model %>% filter(stage_n!="nx") %>% droplevels()
# LUAD_model <- LUAD_model %>% filter(pathologic_stage!="stage iv") %>% droplevels()


#Survival Analysis
#We will use survminer and survival packages in R for this analysis 
#Survival rates for LUSC and LUAD 
KM_LUSC_LUAD <- survfit(Surv(times,patient.vital_status) ~ admin.disease_code,data=LUSC_LUAD_modified)
ggsurvplot(KM_LUSC_LUAD,pval=TRUE, pval.size=3,surv.median.line = "hv",legend="right",title="Survival curves for LUSC and LUAD",legend.title="",legend.labs=c("LUAD","LUSC"),font.main=c(14),font.x=c(12),font.y=c(12),font.tickslab=c(10),xlab="Time in years",risk.table = TRUE,risk.table.y.text=FALSE,risk.table.y.text.col=T,risk.table.height=0.3,risk.table.fontsize=3)
kable(summary(KM_LUSC_LUAD)$table[,7:9],"html",booktabs=T)%>%kable_styling(full_width=TRUE)


# Univariate analysis (LUSC)
### KM survival curves for all the variables 
covariates <- c("gender","agecat","race","smoking_history","pathologic_stage","stage_t","stage_n","stage_m")
KM_LUSC_formula <- sapply(covariates,function(x) as.formula(paste('Surv(time_to_event,vital_status)~',x)))
KM_LUSC_models <- lapply(KM_LUSC_formula,function(x){surv_fit(x,data=LUSC_model)})
KM_LUSC_survdiff <- lapply(KM_LUSC_formula,function(x){survdiff(x,data=LUSC_model)})
KM_LUSC_median <- lapply(KM_LUSC_models,function(x){summary(x)$table[,7:9]})
KM_LUSC_plots <- lapply(KM_LUSC_models,function(x){ggsurvplot(x,pval=TRUE,pval.size=3,surv.median.line = "hv",xlab="Time in Years",font.main=c(14),font.tickslab=c(10),font.x=c(12),font.y=c(12),legend="top",legend.title="",risk.table = TRUE, risk.table.y.text=FALSE, risk.table.col="strata",risk.table.fontsize=3,risk.table.height=0.3,conf.int=TRUE,conf.int.style="step")+guides(colour=guide_legend(nrow=2))})
KM_LUSC_results <- bind_rows(lapply(KM_LUSC_survdiff,function(x){glance(x)}),.id="category")
KM_LUSC_results <- as.data.frame(KM_LUSC_results)
names(KM_LUSC_results)[names(KM_LUSC_results)=="statistic"] <- "chisq" 
kable(KM_LUSC_results,"html",booktabs=T,caption= "Chisq and p-values for LUSC from Kaplan Meier estimates") %>% kable_styling(full_width=TRUE)

#arranging plots 
arrange_ggsurvplots(KM_LUSC_plots[1:2]) 
kable(round(KM_LUSC_median$gender,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
kable(round(KM_LUSC_median$agecat,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="center")

arrange_ggsurvplots(KM_LUSC_plots[3]) 
sh_plot <- ggsurvplot(KM_LUSC_models$smoking_history,pval=TRUE,pval.size=3,surv.median.line = "hv",xlab="Time in Years",font.main=c(14),font.tickslab=c(10),font.x=c(12),font.y=c(12),legend="top",legend.title="",risk.table = TRUE, risk.table.y.text=FALSE, risk.table.col="strata",risk.table.fontsize=3,risk.table.height=0.3,conf.int=TRUE,conf.int.style="step",legend.labs=c("current smoker","reformed smoker","non-smoker"))+guides(colour=guide_legend(nrow=2))
kable(round(KM_LUSC_median$race,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
sh_plot
kable(round(KM_LUSC_median$smoking_history,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="center")

arrange_ggsurvplots(KM_LUSC_plots[5:6])
kable(round(KM_LUSC_median$pathologic_stage,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
kable(round(KM_LUSC_median$stage_t,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
arrange_ggsurvplots(KM_LUSC_plots[7:8])
kable(round(KM_LUSC_median$stage_n,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
kable(round(KM_LUSC_median$stage_m,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="left")

# Cox proportional hazards model (LUSC)
cox_covariates <- c("gender","agecat","race","smoking_history","pathologic_stage","stage_t","stage_n","stage_m")
cox_LUSC_formula <- sapply(cox_covariates,function(x) as.formula(paste('Surv(time_to_event,vital_status)~',x)))
#generating models for all variables 
cox_LUSC_models <- lapply(cox_LUSC_formula,function(x){coxph(x,data=LUSC_model)})
#creating a list of Hazard ratios for variables that have two levels (age,race)
cox_LUSC_results <- lapply(cox_LUSC_models[1:3],function(x){ 
  x <- summary(x)
  p.value<-round(x$wald["pvalue"],digits=2)
  wald.test<- round(x$wald["test"],digits=2)
  beta<-round(x$coef[1],digits=2);#coeficient beta
  HR <-round(x$coef[2],digits = 2);#exp(beta)
  pval <- round(x$coef[5],digits=2);
  HR.confint.lower <- round(x$conf.int[,"lower .95"],digits=2)
  HR.confint.upper <- round(x$conf.int[,"upper .95"],digits=2)
  res<-c(beta, HR, HR.confint.lower,HR.confint.upper,pval,wald.test, p.value)
  names(res)<-c("beta", "HR","95%CI(lower)","95%CI(upper)","p-value","wald.test","wald.p.value")
  return(res)
})

#Testing proportional Hazard assumption
cox_ph <- lapply(cox_LUSC_models,function(x){cox.zph(x,terms=FALSE)})
cox_ph_table <- lapply(cox_ph,function(x){x$table})
cox_ph_table <- do.call(rbind,cox_ph_table)
cox_ph_table <- subset(cox_ph_table,rownames(cox_ph_table) != "GLOBAL")
res_LUSC <- t(as.data.frame(cox_LUSC_results, check.names = FALSE))
res_LUSC <- as.data.frame(res_LUSC)


ggcoxzph(cox_ph$gender)
cox_ph$gender
ggcoxzph(cox_ph$age)
cox_ph$age
ggcoxzph(cox_ph$race)
cox_ph$age
ggcoxzph(cox_ph$pathologic_stage)
cox_ph$pathologic_stage
ggcoxzph(cox_ph$stage_t)
cox_ph$stage_t
ggcoxzph(cox_ph$stage_n)
cox_ph$stage_n
ggcoxzph(cox_ph$stage_m)
cox_ph$stage_m


#creating a list of Hazard ratios and CIs for variables that have multiple levels (cancer stage, staget,n,m,smoking history)
cox_LUSC_sh_ps_HR <- lapply(cox_LUSC_models[4:8],function(x){
  beta <- round(summary(x)$coef[,1],digits=2);
  HR <- round(summary(x)$coef[,2],digits = 2);
  HR.confint.lower <- round(summary(x)$conf.int[,"lower .95"],digits=2)
  HR.confint.upper <- round(summary(x)$conf.int[,"upper .95"],digits=2)
  pval <- round(summary(x)$coef[,5],digits=2);
  res <- cbind(beta,HR,HR.confint.lower,HR.confint.upper,pval)
  #names(res) <- c("beta", "HR","95%CI(lower)","95%CI(upper)","p value")
  return(res)
}
)
#converting the list to a dataframe
#might be a more efficient way to do this but this works for now
cox_LUSC_sh_ps_HR <- as.data.frame(do.call(rbind,cox_LUSC_sh_ps_HR))
rownames(cox_LUSC_sh_ps_HR) <- c("reformed smoker","lifelong nonsmoker","pathologic_stage2","pathologic_stage3","staget2","staget3","staget4","stagen1","stagen2","stagemx")
names(cox_LUSC_sh_ps_HR) <- c("beta","HR","95%CI (L)","95%CI (U)","p value")

#creating a list of p values for variables with multiple levels 
cox_LUSC_sh_ps_pvalues <- lapply(cox_LUSC_models[4:8],function(x){
  p <- round(summary(x)$wald["pvalue"],digits=2)
  t <- round(summary(x)$wald["test"],digits=2)
  res1 <- c(t,p)
  names(res1) <- c("wald.test","p.value")
  return(res1)
})
cox_LUSC_sh_ps_pvalues <- do.call(rbind,cox_LUSC_sh_ps_pvalues)

#printing cox model results as a table 
kable(cox_LUSC_sh_ps_HR,"html",align="l") %>% kable_styling(full_width=FALSE,position="float_left")
kable(cox_LUSC_sh_ps_pvalues,"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="left")

# Univariate Analysis LUAD 
# KM method
KM_LUAD_formula <- sapply(covariates,function(x) as.formula(paste('Surv(time_to_event,vital_status)~',x)))
KM_LUAD_models <- lapply(KM_LUAD_formula,function(x){surv_fit(x,data=LUAD_model)})
KM_LUAD_survdiff <- lapply(KM_LUAD_formula,function(x){survdiff(x,data=LUAD_model)})
KM_LUAD_median <- lapply(KM_LUAD_models,function(x){summary(x)$table[,7:9]})
KM_LUAD_plots <- lapply(KM_LUAD_models,function(x){ggsurvplot(x,pval=TRUE,pval.size=3,surv.median.line = "hv",xlab="Time in Years",font.main=c(14),font.tickslab=c(10),font.x=c(12),font.y=c(12),legend="top",legend.title="", risk.table = TRUE,risk.table.y.text=FALSE, risk.table.col="strata",risk.table.fontsize=3,risk.table.height=0.3,conf.int=TRUE,conf.int.style="step")+guides(colour=guide_legend(nrow=2))})
#KM_LUAD_median <- do.call(rbind,KM_LUAD_median)
#KM_LUAD_median<- as.data.frame(KM_LUAD_median)
#rownames(KM_LUAD_median) <- gsub("^[^=]*=","",rownames(KM_LUAD_median))
KM_LUAD_results <- bind_rows(lapply(KM_LUAD_survdiff,function(x){glance(x)}),.id="category")
KM_LUAD_results <- as.data.frame(KM_LUAD_results)
names(KM_LUAD_results)[names(KM_LUAD_results)=="statistic"] <- "chisq"


#printing KM results as a table 
kable(KM_LUAD_results,"html",booktabs=T) %>%  kable_styling(full_width=TRUE)

#plots and tables for individiual categories that show statistically significan correlation
arrange_ggsurvplots(KM_LUAD_plots[1:2]) 
kable(round(KM_LUAD_median$gender,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="left")
kable(round(KM_LUAD_median$age,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="left")
arrange_ggsurvplots(KM_LUAD_plots[3:4]) 
kable(round(KM_LUAD_median$race,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
kable(round(KM_LUAD_median$smoking_history,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="center")

arrange_ggsurvplots(KM_LUAD_plots[5:6]) 
kable(round(KM_LUAD_median$pathologic_stage,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
kable(round(KM_LUAD_median$stage_t,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="left")

arrange_ggsurvplots(KM_LUAD_plots[7:8]) 
kable(round(KM_LUAD_median$stage_n,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="float_left")
kable(round(KM_LUAD_median$stage_m,digits=2),"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="left")

#calculating pairwise comparisons for statisticallt significant variables that have more than two levels 
pw_LUAD_pathologicstage <- pairwise_survdiff(Surv(time_to_event,vital_status)~pathologic_stage,data=LUAD_model)
symnum(pw_LUAD_pathologicstage$p.value,cutpoints=c(0,0.0001,0.001,0.01,0.05,0.1,1),
       symbols=c("****","***","**","*","+"," "),na="")

pw_LUAD_staget <- pairwise_survdiff(Surv(time_to_event,vital_status)~stage_t,data=LUAD_model)
symnum(pw_LUAD_staget$p.value,cutpoints=c(0,0.0001,0.001,0.01,0.05,0.1,1),
       symbols=c("****","***","**","*","+"," "),na="")

pw_LUAD_stagen <- pairwise_survdiff(Surv(time_to_event,vital_status)~stage_n,data=LUAD_model)
symnum(pw_LUAD_stagen$p.value,cutpoints=c(0,0.0001,0.001,0.01,0.05,0.1,1),
       symbols=c("****","***","**","*","+"," "),na="")


# Cox proportional hazards model LUAD - should ideally write a function for this- maybe at a later time
cox_LUAD_formula <- sapply(cox_covariates,function(x) as.formula(paste('Surv(time_to_event,vital_status)~',x)))
cox_LUAD_models <- lapply(cox_LUAD_formula,function(x){coxph(x,data=LUAD_model)})
cox_LUAD_results <- lapply(cox_LUAD_models[1:3],
                           function(x){ 
                             x <- summary(x)
                             p.value<-round(x$wald["pvalue"],digits=2)
                             wald.test<-round(x$wald["test"],digits=2)
                             beta<-round(x$coef[1],digits=2);#coeficient beta
                             HR <- round(x$coef[2],digits=2);
                             pval <- round(x$coef[,5],digits=2);
                             HR.confint.lower <- round(x$conf.int[,"lower .95"],digits=2)
                             HR.confint.upper <- round(x$conf.int[,"upper .95"],digits=2)
                             res<-c(beta, HR, HR.confint.lower,HR.confint.upper,pval,wald.test, p.value)
                             names(res)<-c("beta", "HR","95%_CI_lower","95%_CI_upper","p value","wald.test", "wald.p.value")
                             return(res)
                           })
#Testing for Cox proportional hazard assumption
cox_ph_LUAD <- lapply(cox_LUAD_models,function(x){cox.zph(x,terms=FALSE)})
cox_ph_LUAD_table <- lapply(cox_ph_LUAD,function(x){x$table})
cox_ph_LUAD_table <- do.call(rbind,cox_ph_LUAD_table)
cox_ph_LUAD_table <- subset(cox_ph_LUAD_table,rownames(cox_ph_LUAD_table) != "GLOBAL")
res_LUAD <- t(as.data.frame(cox_LUAD_results, check.names = FALSE))
res_LUAD <- as.data.frame(res_LUAD)

ggcoxzph(cox_ph_LUAD$gender)
cox_ph_LUAD$gender
ggcoxzph(cox_ph_LUAD$age)
cox_ph_LUAD$age
ggcoxzph(cox_ph_LUAD$race)
cox_ph_LUAD$age
ggcoxzph(cox_ph_LUAD$pathologic_stage)
cox_ph_LUAD$pathologic_stage
ggcoxzph(cox_ph_LUAD$stage_t)
cox_ph_LUAD$stage_t
ggcoxzph(cox_ph_LUAD$stage_n)
cox_ph_LUAD$stage_n
ggcoxzph(cox_ph_LUAD$stage_m)
cox_ph_LUAD$stage_m

#printing results from proportional hazard assumption tests
kable(res_LUAD,"html",booktabs=T,caption="Univariate analysis using Cox proportion hazards") %>%  kable_styling(full_width=TRUE)

cox_LUAD_sh_ps_HR <- lapply(cox_LUAD_models[4:8],function(x){
  beta <- round(summary(x)$coef[,1],digits=2);
  HR <- round(summary(x)$coef[,2],digits = 2);
  pval <- round(summary(x)$coef[,5],digits=2);
  HR.confint.lower <- round(summary(x)$conf.int[,"lower .95"],digits=2)
  HR.confint.upper <- round(summary(x)$conf.int[,"upper .95"],digits=2)
  res <- cbind(beta,HR,HR.confint.lower,HR.confint.upper,pval)
  #names(res) <- c("beta", "HR","95%_CI_lower","95%_CI_upper","p-value")
  return(res)
}
)
cox_LUAD_sh_ps_HR <- as.data.frame(do.call(rbind,cox_LUAD_sh_ps_HR))
rownames(cox_LUAD_sh_ps_HR) <- c("reformed_smoker","non-smoker","pathologic_stage2","pathologic_stage3","pathologic_stage4","staget2","staget3","staget4","stagen1","stagen2","stagem1","stagemx")
names(cox_LUAD_sh_ps_HR) <- c("beta", "HR","0.95LCL","0.95UCL","P-value")

cox_LUAD_sh_ps_pvalues <- lapply(cox_LUAD_models[4:8],function(x){
  p <- round(summary(x)$wald["pvalue"],digits = 2)
  t <- round(summary(x)$wald["test"],digits=2)
  res1 <- c(t,p)
  names(res1) <- c("wald.test","wald.p.value")
  return(res1)
})
cox_LUAD_sh_ps_pvalues <- do.call(rbind,cox_LUAD_sh_ps_pvalues)
#printing results from cox models for lUAD
kable(cox_LUAD_sh_ps_HR,"html",align="l") %>% kable_styling(full_width=FALSE,position="float_left") %>% row_spec(0,angle=0)
kable(cox_LUAD_sh_ps_pvalues,"html",booktabs=T) %>% kable_styling(full_width=FALSE,position="left")


#Multivariate analysis for LUAD 
# Model 1 (agecat,gender,pathologic_stage,smoking history)
# Testing the proportional hazards assumption

cox_multivariate1 <- coxph(Surv(time_to_event,vital_status)~agecat+gender+pathologic_stage+smoking_history,data=LUAD_model)
model1_zph <- cox.zph(cox_multivariate1)
#table and forest plot for model 1
ggcoxzph(model1_zph)
ggforest(cox_multivariate1,data=LUAD_model)

  
#Model 2 - Includes age, tnm stages t,n and smoking history
cox_multivariate2 <- coxph(Surv(time_to_event,vital_status)~agecat+gender+stage_t+stage_n+stage_m+smoking_history,data=LUAD_model)
model2 <- cox.zph(cox_multivariate2)
#table and forest plot for model2
ggcoxzph(model2)
ggforest(cox_multivariate2,data=LUAD_model)
