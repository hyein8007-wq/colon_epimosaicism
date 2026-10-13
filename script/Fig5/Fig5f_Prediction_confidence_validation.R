################################################################################
# Figure 5f. Validation of phylogeny refinement prediction confidence
################################################################################

# 1. Libraries -----------------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(scales)

# 2. Input data ----------------------------------------------------------------
validation <- read_tsv("input/Fig5/phylogeny_refinement_validation.tsv",show_col_types=FALSE)
stopifnot(all(c("Shared somatic mutations","Top support","Second support","Confidence","Outcome") %in% names(validation)))

# 3. Validation data processing ------------------------------------------------
# Confidence is recalculated from unrounded support counts rather than the rounded table value.
df2 <- validation %>% filter(!is.na(`Shared somatic mutations`),`Shared somatic mutations`>5,Outcome %in% c("Correct","Incorrect")) %>% mutate(confidence=`Top support`/(`Top support`+`Second support`),correct_num=as.numeric(Outcome=="Correct")) %>% filter(is.finite(confidence))
stopifnot(nrow(df2)>0,length(unique(df2$correct_num))==2,length(unique(df2$confidence))>1)

# 4. Logistic regression -------------------------------------------------------
fit <- glm(correct_num~confidence,family=binomial(),data=df2)
pval <- summary(fit)$coefficients["confidence","Pr(>|z|)"]
newdat <- data.frame(confidence=seq(min(df2$confidence),max(df2$confidence),length.out=300))
pred <- predict(fit,newdata=newdat,type="link",se.fit=TRUE)
newdat <- newdat %>% mutate(fit=plogis(pred$fit),lower=plogis(pred$fit-1.96*pred$se.fit),upper=plogis(pred$fit+1.96*pred$se.fit))

# 5. Figure 5f -----------------------------------------------------------------
df_plot <- df2 %>% mutate(conf_round=round(confidence,2)) %>% group_by(correct_num,conf_round) %>% mutate(stack=row_number()) %>% ungroup() %>% mutate(y_raw=ifelse(correct_num==1,1.055+(stack-1)*0.025,-0.055-(stack-1)*0.025))
p <- ggplot()+
  geom_ribbon(data=newdat,aes(confidence,ymin=lower,ymax=upper),fill="grey70",alpha=0.4)+
  geom_line(data=newdat,aes(confidence,fit),linewidth=1.2)+
  geom_point(data=df_plot,aes(confidence,y_raw,color=factor(correct_num)),size=2.1,show.legend=FALSE)+
  scale_color_manual(values=c("0"="#E88D95","1"="#5899D4"))+
  annotate("text",x=0.505,y=1.115,label="Correct",hjust=0,size=3.7)+
  annotate("text",x=0.505,y=-0.115,label="Incorrect",hjust=0,size=3.7)+
  annotate("text",x=0.73,y=0.22,hjust=0,size=3.8,label=paste0("Logistic regression\nP = ",formatC(pval,format="f",digits=3)))+
  scale_x_continuous(limits=c(0.49,1.02),breaks=seq(0.5,1,0.1),expand=expansion(mult=0))+
  scale_y_continuous(limits=c(-0.16,1.16),breaks=seq(0,1,0.2),labels=percent_format(accuracy=1),expand=expansion(mult=0))+
  labs(x="Prediction confidence",y="Probability of correct prediction")+
  coord_cartesian(clip="off")+theme_classic(base_size=14)+
  theme(axis.line=element_line(linewidth=0.8),axis.ticks=element_line(linewidth=0.8),axis.text=element_text(color="black"),plot.margin=margin(15,15,15,15))
print(p)
