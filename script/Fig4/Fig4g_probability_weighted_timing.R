################################################################################
# Figure 4g. Probability-weighted inference of XCI timing
# Run from the repository root. No output files are written.
################################################################################

# 1. Libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# 2. Input data
input_dir <- "input/Fig4"
hc50_raw <- read_tsv(file.path(input_dir,"HC50_RandomSample_100_final_res_260919_v20.tsv"),show_col_types=FALSE)
hcs_raw <- read_tsv(file.path(input_dir,"HCs_final_res_260919_v20.tsv"),show_col_types=FALSE)
stopifnot(all(c("n_subtree","mutation_point","weight") %in% names(hc50_raw)),all(c("sample_id","mutation_point","weight") %in% names(hcs_raw)))

# 3. HC50: normalize each random subtree, then average as one donor
hc50_subtree <- hc50_raw %>% group_by(n_subtree,mutation_point) %>% summarise(weight=sum(weight,na.rm=TRUE),.groups="drop") %>% group_by(n_subtree) %>% mutate(weight_norm=weight/sum(weight,na.rm=TRUE)) %>% ungroup()
hc50_subtree_complete <- hc50_subtree %>% select(n_subtree,mutation_point,weight_norm) %>% complete(n_subtree,mutation_point=seq(min(mutation_point,na.rm=TRUE),max(mutation_point,na.rm=TRUE)),fill=list(weight_norm=0))
hc50_donor <- hc50_subtree_complete %>% group_by(mutation_point) %>% summarise(weight=mean(weight_norm),.groups="drop") %>% mutate(sample_id="HC50")

# 4. Other donors: normalize each donor independently
hcs_donor <- hcs_raw %>% group_by(sample_id,mutation_point) %>% summarise(weight=sum(weight,na.rm=TRUE),.groups="drop") %>% group_by(sample_id) %>% mutate(weight=weight/sum(weight,na.rm=TRUE)) %>% ungroup()
stopifnot(!"HC50" %in% hcs_donor$sample_id)

# 5. Equal-weight donor-level probability distribution
mutation_density <- bind_rows(hcs_donor,hc50_donor)
mutation_density_complete <- mutation_density %>% complete(sample_id,mutation_point=seq(min(mutation_point,na.rm=TRUE),max(mutation_point,na.rm=TRUE)),fill=list(weight=0))
mutation_density_total <- mutation_density_complete %>% group_by(mutation_point) %>% summarise(weight=mean(weight),.groups="drop")
probability_check <- mutation_density %>% group_by(sample_id) %>% summarise(total_probability=sum(weight),.groups="drop")
print(probability_check)
stopifnot(all(abs(probability_check$total_probability-1)<1e-6),abs(sum(mutation_density_total$weight)-1)<1e-6)

# 6. Plot: donor KDEs, pooled KDE and density strip
xmin <- -7
xmax <- 30
g <- ggplot() +
  geom_density(data=mutation_density_complete,aes(x=mutation_point,weight=weight,group=sample_id),color="gray",alpha=0.5,linewidth=1,adjust=1,from=xmin,to=xmax) +
  geom_density(data=mutation_density_total,aes(x=mutation_point,weight=weight),color="blue",linewidth=1.5,adjust=1,from=xmin,to=xmax) +
  geom_tile(data=mutation_density_total,aes(x=mutation_point,y=-0.003,weight=weight,fill=after_stat(density)),stat="density",height=0.006,width=0.15,adjust=1,from=xmin,to=xmax) +
  scale_fill_viridis_c(option="C") +
  scale_x_continuous(limits=c(xmin,xmax),breaks=seq(-5,30,5),expand=c(0,0)) +
  labs(x="Number of somatic mutations",y="Probability density",fill="Density") + theme_classic()
print(g)

# 7. KDE peak and highest-density intervals (original calculation)
d <- density(mutation_density_total$mutation_point,weights=mutation_density_total$weight/sum(mutation_density_total$weight),adjust=1,from=xmin,to=xmax,n=4096)
peak_x <- d$x[which.max(d$y)]
peak_density <- max(d$y)
dx <- d$x[2]-d$x[1]
kde_df <- tibble(x=d$x,density=d$y,mass=d$y*dx)
tmp <- kde_df %>% arrange(desc(density)) %>% mutate(cum_mass=cumsum(mass))
cutoff_50 <- tmp %>% filter(cum_mass>=0.50) %>% slice(1) %>% pull(density)
hdi50 <- kde_df %>% filter(density>=cutoff_50)
hdi50_range <- c(lower=min(hdi50$x),upper=max(hdi50$x))
cutoff_95 <- tmp %>% filter(cum_mass>=0.95) %>% slice(1) %>% pull(density)
hdi95 <- kde_df %>% filter(density>=cutoff_95)
hdi95_range <- c(lower=min(hdi95$x),upper=max(hdi95$x))
cat("Peak:",round(peak_x,2),"\n","50% HDI:",round(hdi50_range[1],2),"-",round(hdi50_range[2],2),"\n","95% HDI:",round(hdi95_range[1],2),"-",round(hdi95_range[2],2),"\n")
