library(tidyverse)
library(haven)
library(rspiro)

setwd("~/Dropbox/My Mac (Emily’s MacBook Air)/Documents/Fellowship/MSCE/DataScience/BMIN503_Final_Project")

set.seed(36)

# Load spirometry
spirometry1 <- read_sas("data/spxraw_e.sas7bdat")
spirometry2 <- read_sas("data/spxraw_f.sas7bdat")
spirometry3 <- read_sas("data/spxraw_g.sas7bdat")

all_spirometry <- bind_rows(spirometry1, spirometry2, spirometry3)
rm(spirometry1,spirometry2,spirometry3)

all_spirometry <- all_spirometry %>% rename(seqn = SEQN)

# Clean spirometry

all_spirometry <- 
  all_spirometry %>%
  mutate(length = str_count(SPXRAW, "\\d,"))

all_spirometry %>%
  ggplot(aes(x=length)) + geom_histogram(binwidth = 1)

# Load body measurements
bm1 <- read_xpt("data/BMX_E.XPT") %>% select(SEQN, BMXWT, BMXHT)
bm2 <- read_xpt("data/BMX_F.XPT") %>% select(SEQN, BMXWT, BMXHT)
bm3 <- read_xpt("data/BMX_G.XPT") %>% select(SEQN, BMXWT, BMXHT)

all_bm <- bind_rows(bm1, bm2, bm3) %>% rename(seqn = SEQN)
rm(bm1, bm2, bm3)

all_bm <- all_bm %>%
  rename(wt = BMXWT,
         ht = BMXHT)


# Load mortality follow up data
# read in the fixed-width format ASCII file
mort1 <- read_fwf(file="data/NHANES_2007_2008_MORT_2019_PUBLIC.dat",
                col_types = "iiiiiiii",
                fwf_cols(seqn = c(1,6),
                         eligstat = c(15,15),
                         mortstat = c(16,16),
                         ucod_leading = c(17,19),
                         diabetes = c(20,20),
                         hyperten = c(21,21),
                         permth_int = c(43,45),
                         permth_exm = c(46,48)
                ),
                na = c("", ".")
)

mort2 <- read_fwf(file="data/NHANES_2009_2010_MORT_2019_PUBLIC.dat",
                  col_types = "iiiiiiii",
                  fwf_cols(seqn = c(1,6),
                           eligstat = c(15,15),
                           mortstat = c(16,16),
                           ucod_leading = c(17,19),
                           diabetes = c(20,20),
                           hyperten = c(21,21),
                           permth_int = c(43,45),
                           permth_exm = c(46,48)
                  ),
                  na = c("", ".")
)

mort3 <- read_fwf(file="data/NHANES_2011_2012_MORT_2019_PUBLIC.dat",
                  col_types = "iiiiiiii",
                  fwf_cols(seqn = c(1,6),
                           eligstat = c(15,15),
                           mortstat = c(16,16),
                           ucod_leading = c(17,19),
                           diabetes = c(20,20),
                           hyperten = c(21,21),
                           permth_int = c(43,45),
                           permth_exm = c(46,48)
                  ),
                  na = c("", ".")
)

all_mortality <- bind_rows(mort1, mort2, mort3)
rm(mort1, mort2, mort3)

# Load demographic data
demo1 <- read_xpt("data/DEMO_E.XPT")
demo2 <- read_xpt("data/DEMO_F.XPT")
demo3 <- read_xpt("data/DEMO_G.XPT")

# The demographic variables are different across NHANES years
# so for now I will create a subset with only variables of interest to me

demo1 <- demo1 %>% select(SEQN, RIDAGEEX, RIDRETH1, RIAGENDR)
demo2 <- demo2 %>% select(SEQN, RIDAGEEX, RIDRETH1, RIAGENDR)
demo3 <- demo3 %>% select(SEQN, RIDEXAGM, RIDRETH1, RIAGENDR) %>%
  rename(RIDAGEEX = RIDEXAGM)

all_demo <- bind_rows(demo1, demo2, demo3)
rm(demo1, demo2, demo3)

all_demo <-
  all_demo %>% rename(
  seqn = SEQN,
  age_months = RIDAGEEX,
  gender = RIAGENDR,
  race = RIDRETH1
)

all_demo <- all_demo %>%
  mutate(age_years = age_months/12)

# Load PFT measurements
pft1 <- read_xpt("data/SPX_E.xpt")
pft2 <- read_xpt("data/SPX_F.xpt")
pft3 <- read_xpt("data/SPX_G.xpt")

all_pft <- bind_rows(pft1, pft2, pft3)
rm(pft1, pft2, pft3)

all_pft <- all_pft %>% rename(seqn = SEQN)


# Now to perform the laborious task of generating
# predicted spirometry values for these people
# so that I can measure BD responsiveness

# pred_GLI takes age in years, height in m,
# and gender as 1=male 2=female (c/w NHANES)

bd_response <- 
  all_pft %>% left_join(all_demo, by = "seqn") %>%
  left_join(all_bm, by = "seqn") %>%
  filter(SPXBSTAT==1) %>% # Filter on those who performed BD testing
  filter(age_years >= 18) %>%
  select(seqn, SPXNFEV1, SPXNFVC, SPXBFEV1, SPXBFVC ,age_years, gender, race, ht) %>%
  mutate(ht = ht/100) %>% # Convert height to meters
  rename(fev1 = SPXNFEV1) %>%
  mutate(fev1 = fev1/1000) %>% # Convert to L
  rename(fvc = SPXNFVC) %>%
  mutate(fvc = fvc/1000) %>% # Convert to L
  rename(fev1_post = SPXBFEV1) %>%
  mutate(fev1_post = fev1_post/1000) %>%
  rename(fvc_post = SPXBFVC) %>%
  mutate(fvc_post = fvc_post/1000) %>%
  mutate(fev1_pred =
           pred_GLIgl(age_years,
                    ht,
                    gender,
                    param = "FEV1")) %>%
  mutate(fvc_pred =
           pred_GLIgl(age_years,
                    ht,
                    gender,
                    param = "FVC")) %>%
  mutate(fev1_change = (fev1_post-fev1)/fev1_pred,
         fvc_change = (fvc_post-fvc)/fvc_pred,
         bd_response = case_when(
           fev1_change >= 0.1 ~ TRUE,
           fvc_change >= 0.1 ~ TRUE,
           .default = FALSE
         )
  )

bd_response %>% count(bd_response)





all_spirometry_cleaned <-
  all_spirometry %>% left_join(all_demo, by = "seqn") %>%
  filter(age_years >= 18) %>%
  filter(SPATTYPE=="Pre") %>%
  filter(SPAPLAT=="Y") %>%
  filter(SPAACC=="Y") %>%
  filter(SPAQEFF=="A" | SPAQEFF=="B") %>%
  left_join(all_mortality, by="seqn") %>%
  select(seqn, mortstat,age_years,race,gender,SPXRAW,SPAMANU) %>%
  filter(!is.na(mortstat))



regression_data <- 
  all_pft %>% left_join(all_demo, by = "seqn") %>%
  left_join(all_mortality, by = "seqn") %>%
  select(seqn, gender, age_months, SPXNFEV1, SPXNFVC, SPXNF257, race,
         mortstat) %>%
  rename(fev1 = SPXNFEV1,
         fvc = SPXNFVC,
         fef257 = SPXNF257) %>%
  mutate(race = factor(race)) %>%
  mutate(gender = case_when(
    gender==1 ~ "male",
    gender==2 ~ "female"
  )) %>%
  mutate(gender = factor(gender)) %>%
  mutate(mortstat = factor(mortstat)) %>%
  mutate(assignment = rbinom(nrow(regression_data),1,0.8))

train_data <- filter(regression_data, assignment==1)
val_data <- filter(regression_data, assignment==0)

m0 <- glm(gender ~ fev1 + fvc + fef257 + age_months + race,
          data = train_data,
          family = "binomial")

summary(m0)

predictions <- 
  predict.glm(
    m0,
    newdata = val_data,
    type = "response"
  )

predictions <-
  bind_cols(predictions, val_data) %>%
  rename(prediction = `...1`) %>%
  mutate(prediction = case_when(
    prediction > 0.5 ~ "female",
    prediction <= 0.5 ~ "male"
  )) %>%
  filter(!is.na(gender) & !is.na(prediction)) %>%
  mutate(accurate = gender==prediction)

predictions %>%
  summarise(sum(accurate)/n())

table(predictions$prediction, predictions$gender)

# Export for loading into Python for ML
spiro_for_py <-
  all_spirometry_cleaned %>%
  select(seqn, SPXRAW, gender)

spiro_for_py %>% mutate(length = str_count(SPXRAW, "\\d,")) %>% arrange(length) %>%
  ggplot(aes(x=length)) + geom_histogram(binwidth=1)

write_csv(spiro_for_py, "data/spiro_for_py.csv")

# Split the SPXRAW strings into lists of numeric values
all_spirometry_cleaned$SPXRAW_numeric <- 
  lapply(all_spirometry_cleaned$SPXRAW, function(x) as.numeric(unlist(strsplit(x, ","))))

plot_data <- do.call(rbind, lapply(1:nrow(all_spirometry_cleaned), function(i) {
  data.frame(
    x = 1:length(all_spirometry_cleaned$SPXRAW_numeric[[i]]),
    y = all_spirometry_cleaned$SPXRAW_numeric[[i]],
    seqn = all_spirometry_cleaned$seqn[i],  # Retain the seqn variable for each row
    gender = all_spirometry_cleaned$gender[i], # Retain the gender variable
    trial = all_spirometry_cleaned$SPAMANU[i] # Retain the trial variable
  )
}))

plot_data %>%
  mutate(gender = factor(gender)) %>%
  filter(x %% 10 == 0) %>%
  group_by(gender,x) %>%
  summarise(
    y_mean = quantile(y,99/100)
  ) %>%
  ggplot(aes(x=x, y=y_mean, color=gender)) +
  geom_line()

plot_data %>%
  filter(x%%100==0) %>%
  mutate(gender = factor(gender)) %>%
  ggplot(aes(x=x, y=y, color=gender)) +
  geom_point(alpha=0.1) +
  theme_bw()

#MORTSTAT: Final Mortality Status
# 0 = Assumed alive
# 1 = Assumed deceased
# <NA> = Ineligible or under age 18

#UCOD_LEADING: Underlying Cause of Death: Recode
# 1 = Diseases of heart (I00-I09, I11, I13, I20-I51)
# 2 = Malignant neoplasms (C00-C97)
# 3 = Chronic lower respiratory diseases (J40-J47)
# 4 = Accidents (unintentional injuries) (V01-X59, Y85-Y86)
# 5 = Cerebrovascular diseases (I60-I69)
# 6 = Alzheimer's disease (G30)
# 7 = Diabetes mellitus (E10-E14)
# 8 = Influenza and pneumonia (J09-J18)
# 9 = Nephritis, nephrotic syndrome and nephrosis (N00-N07, N17-N19, N25-N27)
# 10 = All other causes (residual)
# <NA> = Ineligible, under age 18, assumed alive, or no cause of death data available