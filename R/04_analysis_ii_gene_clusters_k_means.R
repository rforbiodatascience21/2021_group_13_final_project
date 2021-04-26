#I know the below data wrangling is ugly

gravier_data <- load("data/_raw/gravier.RData")

gravier_data <-
  bind_cols(gravier %>% 
              pluck("y") %>% 
              as.tibble(),gravier %>% 
              pluck("x") %>% 
              as.tibble()) %>% 
  as.tibble() %>% 
  mutate(outcome=case_when(value=="good"~0,value=="poor"~1)) %>% 
  select(outcome,everything(),-value) 
model1 <- gravier_data %>% glm(outcome~g2E09,data=.,family=binomial(link="logit"))

gravier_data_long <- 
  gravier_data %>% 
  pivot_longer(-outcome,names_to="gene",values_to = "log2_expr_level") %>% 
  mutate(log2_expr_level=round(log2_expr_level,5))

gravier_data_nested <- gravier_data_long %>% group_by(gene) %>% nest() %>% ungroup()
rm(model1,gravier_data_long)

set.seed(42)
selected <- sample_n(gravier_data_nested,100)
gravier_data_nested_long <- 
  selected %>% 
  mutate(mdl=map(data,~glm(outcome~log2_expr_level,data=.x,family=binomial(link="logit"))))
rm(gravier_data_nested,selected)
gravier_data_nested_long <-
  gravier_data_nested_long %>% 
  mutate(tidy=map(mdl,tidy,conf.int=TRUE)) %>% 
  unnest(tidy) %>% 
  filter(term!="(Intercept)") %>% 
  mutate(identified_as=case_when(p.value<0.05~"significant",TRUE~"unsignificant"),neg_log_p=-log10(p.value))

#select outcome and gene expression levels for each patient
gravier_data_wide = gravier_data %>%
  select(outcome, pull(gravier_data_nested_long, gene))


#mutate to get better readable labels
gravier_data_wide <-
  gravier_data_wide %>% 
  mutate(outcome=as.factor(outcome)) %>% 
  mutate(outcome=case_when(outcome=="0"~"good",TRUE~"poor"))

#k-means gravier-wide: maybe look for clusters of expression levels for good and bad outcome, first with all genes, then maybe just some PCs

clusters <- 
  gravier_data_wide %>% 
  select(-outcome) %>% 
  kmeans(centers=2)

#way to check if simple 2 centres clustering can predict outcome:
augment(clusters,gravier_data_wide) %>% select(outcome,.cluster) %>% mutate(pred=case_when((outcome=="good" & .cluster==1) | (outcome=="poor" & .cluster==2) ~ "correct",TRUE~"wrong")) %>% count(pred=="correct")

#--> 81 incorrectly, 87 corretly predicted

#now try different numbers of clusters, maybe they can predict the outcome better as there may exist more than 2 different genotypes with the same outcome
data_to_cluster <-
  gravier_data_wide %>% 
  select(-outcome)
kclusters <- 
  tibble(k=1:100) %>% 
  mutate(cluster=map(k,~kmeans(data_to_cluster,.x)),
         tidy=map(cluster,tidy),
         glance=map(cluster,glance),
         augmented=map(cluster,augment,gravier_data_wide))

#take the augmented data and try to see if a certain amount of clusters is better
#idea: for each number of centres calculate the ratio of good/poor in each cluster, add those ratios up and divide by the number of clusters to get an average of how good the ratio of classification was for each model(with x amount of clusters)
kclusters %>% 
  select(k,augmented) %>% 
  unnest(cols=c(augmented)) %>% 
  select(k,outcome,.cluster) %>% 
  group_by(k,.cluster,outcome) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = outcome,values_from = n) %>%
  ungroup() %>% 
  mutate(poor=case_when(is.na(poor)~as.integer(0),TRUE~poor),good=case_when(is.na(good)~as.integer(0),TRUE~good)) %>% 
  mutate(ratio=case_when(poor>good~good/poor,TRUE~poor/good)) %>% 
  group_by(k) %>% 
  summarise(sum=sum(ratio/k)) %>% 
  arrange(sum)

# Clean up ----------------------------------------------------------------

rm(clusters,data_to_cluster,kclusters,gravier,gravier_data,gravier_data_wide,gravier_data_nested_long)
