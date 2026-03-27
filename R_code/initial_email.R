library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyverse)
library(skimr)
library(inspectdf)
library(table1)
library(naniar)
library(knitr)

dataset <-
  read_excel(here::here("raw-data", "FITFR-final-data-v3.xlsx"),
             skip = 1,
             n_max = 100)

head(dataset)

summary(dataset$age, useNA = "ifany")
table(dataset$sex, useNA = "ifany")
table(dataset$treatment, useNA = "ifany")
table(dataset$smoker, useNA = "ifany")
summary(dataset$AF)
table(dataset$PR)

snp_cols <- grep("^rs", names(dataset), value = TRUE)

# valid_genotypes <- c(
#   "AA","AC","AG","AT",
#   "CC","CG","CT",
#   "GG","GT",
#   "TT"
# )
# invalid_counts <- sapply(dataset[snp_cols],
#       function(col) {sum(!col %in% valid_genotypes & !is.na(col))})
# invalid_counts

dataset[snp_cols] <- lapply(dataset[snp_cols], as.factor)

dataset <- dataset %>%
  mutate(
    age = as.integer(age),
    age = ifelse(age >= 0 & age <= 120, age, NA),
    sex = ifelse(tolower(sex) %in% c("male", "m"),
                 "Male",
                 ifelse(tolower(sex) %in% c("female", "f"), "Female", NA)),
    treatment = toupper(treatment),
    smoker = ifelse(tolower(smoker) %in% c("n"), "N",
                    ifelse(tolower(smoker) %in% c("y"), "Y", NA))
    )

miss_var_summary(dataset)


dataset <- dataset %>%
  filter(complete.cases(.))

dataset$sex <- as.factor(dataset$sex)
dataset$treatment <- as.factor(dataset$treatment)
dataset$smoker <- as.factor(dataset$smoker)
dataset$PR <- as.factor(dataset$PR)


summary(dataset$age, useNA = "ifany")
table(dataset$sex, useNA = "ifany")
table(dataset$treatment, useNA = "ifany")
table(dataset$smoker, useNA = "ifany")
summary(dataset$AF)
table(dataset$PR)

to_remove <- dataset %>%
  select(where(~ length(unique(.)) < 2))

dataset <- dataset %>%
  select(where(~ length(unique(.)) >= 2))

dataset$PR_num <- as.numeric(as.factor(dataset$PR)) - 1

saveRDS(dataset, here::here("data", "initial_dataset_clean.rds"))

dataset <- readRDS(here::here("data", "initial_dataset_clean.rds"))

snp_cols <- grep("^rs", names(dataset), value = TRUE)

# Plots and summary statistics
############################################################

AF_summary_stats_groupedBy_treatment <- dataset %>%
  group_by(treatment) %>%
  summarise(
      mean = round(mean(AF),2),
      sd   = round(sd(AF),2),
      min  = round(min(AF),2),
      max  = round(max(AF),2)
  )
kable(AF_summary_stats_groupedBy_treatment,
      caption = "AF Summary by Treatment")

prop.table(table(dataset$treatment, dataset$PR))

#Prop table grouped by treatment
prop.table(table(dataset$treatment, dataset$PR), margin = 1)


# AF vs Age by Treatment [No trend line]
ggplot(dataset, aes(x = age, y = AF, colour = treatment)) +
  geom_point() +
  labs(title = "AF vs Age by Treatment",
       x = "Age",
       y = "AF")

# AF vs Age by Treatment [Trend line]
ggplot(dataset, aes(x = age, y = AF, colour = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = False) +
  labs(title = "AF vs Age by Treatment",
       x = "Age",
       y = "AF")

# AF vs Age by Treatment [Separate plots]
ggplot(dataset, aes(x = age, y = AF)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~treatment) +
  labs(title = "AF vs Age by Treatment",
       x = "Age",
       y = "AF")


# AF vs Treatment
ggplot(dataset, aes(x = treatment, y = AF, fill = treatment)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "AF by Treatment",
       x = "Treatment",
       y = "AF")


# Bar counts PR by treatment
ggplot(dataset, aes(x = treatment, fill = PR)) +
  geom_bar() +
  labs(title = "PR distribution by Treatment",
       x = "Treatment",
       y = "Count")

# Bar counts PR by sex
ggplot(dataset, aes(x = sex, fill = PR)) +
  geom_bar() +
  labs(title = "PR distribution by Sex",
       x = "Sex",
       y = "Count")

# Bar counts PR by treatment (counts in figure)
ggplot(dataset, aes(x = treatment, fill = PR)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_stack(vjust = 0.5)) +
  labs(title = "PR distribution by Treatment",
       x = "Treatment",
       y = "Count")


# Bar proportions PR by treatment
ggplot(dataset, aes(x = treatment, fill = PR)) +
  geom_bar(position = "fill") +
  labs(title = "PR Proportion by Treatment",
       x = "Treatment",
       y = "Proportion")


# PR counts by treatment level
PR_counts_by_treatment_level <- dataset %>%
  group_by(treatment, PR) %>%
  summarise(n = n(), .groups = "drop")
PR_counts_by_treatment_level


# AF vs PR
ggplot(dataset, aes(x = PR, y = AF, fill = PR)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "AF by PR",
       x = "PR",
       y = "AF")


# PR vs age
dataset$PR_num <- as.numeric(as.factor(dataset$PR)) - 1
ggplot(dataset, aes(x = age, y = PR_num)) +
  geom_jitter(height = 0.05) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"), se = FALSE) +
  scale_y_continuous(
    breaks = c(0,1),
    labels = levels(dataset$PR)
  ) +
  labs(title = "PR vs Age",
       x = "Age",
       y = "PR")



#SNPs ###
######### PR

# Bar counts PR by rs91154 (counts in figure)
ggplot(dataset, aes(x = rs91154, fill = PR)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_stack(vjust = 0.5)) +
  labs(title = "PR distribution by rs91154",
       x = "rs91154",
       y = "Count")

# Bar proportions PR by rs91154
ggplot(dataset, aes(x = rs91154, fill = PR)) +
  geom_bar(position = "fill") +
  labs(title = "PR Proportion by rs91154",
       x = "rs91154",
       y = "Proportion")


# Bar counts PR by rs56579 (counts in figure)
ggplot(dataset, aes(x = rs56579, fill = PR)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_stack(vjust = 0.5)) +
  labs(title = "PR distribution by rs56579",
       x = "rsrs56579",
       y = "Count")

# Bar proportions PR by rs56579
ggplot(dataset, aes(x = rs56579, fill = PR)) +
  geom_bar(position = "fill") +
  labs(title = "PR Proportion by rs56579",
       x = "rsrs56579",
       y = "Proportion")



######### AF
ggplot(dataset, aes(x = rs57179, y = AF, fill = rs57179)) +
  geom_boxplot() +
  labs(title = "AF by rs57179",
       x = "rs57179",
       y = "AF")


ggplot(dataset, aes(x = rs19486, y = AF, fill = rs19486)) +
  geom_boxplot() +
  labs(title = "AF by rs19486",
       x = "rs19486",
       y = "AF")


# Models
##########################################################
# AF_formula <- as.formula(
#   paste("AF ~ age + treatment + sex +",
#         paste(snp_cols, collapse = " + "))
# )
#
# PR_formula <- as.formula(
#   paste("PR ~ age + treatment + sex +",
#         paste(snp_cols, collapse = " + "))
# )


fit_bivariate <- function(x, y, df) {
  df <- df |>
    select(all_of(c(x, y))) |>
    na.omit()
  model_form <- reformulate(x, y)
  null_form <- formula(str_glue("{y} ~ 1"))
  full <- lm(model_form, df)
  null <- lm(null_form, df)
  pv <- anova(full, null, test = "LRT") |>
    broom::tidy() |>
    slice(2) |>
    pull(p.value)
  tibble(term = x, pv = pv)
}
fit_bivariates <- function(y, preds, df) {
  bivariates <-
    preds |>
    map(safely(\(x) fit_bivariate(x, y, df))) |>
    transpose()
  bivariates$result |>
    list_rbind() |>
    mutate(adj_p = p.adjust(pv, method = "fdr")) |>
    filter(!is.na(adj_p)) |>
    mutate(outcome = y) |>
    arrange(adj_p)
}

preds <- c("age", "sex", "treatment", "PR", snp_cols)
AF_model <- fit_bivariates("AF", preds, dataset)


fit_bivariateLogistic <- function(x, y, df) {
  df <- df |>
    select(all_of(c(x, y))) |>
    na.omit()
  model_form <- reformulate(x, y)
  null_form <- formula(str_glue("{y} ~ 1"))
  full  <- glm(model_form, data = df, family = binomial)
  null  <- glm(null_form, data = df, family = binomial)
  pv <- anova(full, null, test = "LRT") |>
    broom::tidy() |>
    slice(2) |>
    pull(p.value)
  tibble(term = x, pv = pv)
}
fit_bivariatesLogistic <- function(y, preds, df) {
  bivariates <-
    preds |>
    map(safely(\(x) fit_bivariateLogistic(x, y, df))) |>
    transpose()
  bivariates$result |>
    list_rbind() |>
    mutate(adj_p = p.adjust(pv, method = "fdr")) |>
    filter(!is.na(adj_p)) |>
    mutate(outcome = y)|>
    arrange(adj_p)
}

preds_logistic <- c("age", "sex", "treatment", "AF", snp_cols)
PR_model <- fit_bivariatesLogistic("PR", preds_logistic, dataset)


# Other inspection
#############################################
skim(dataset)

# Numeric columns
inspect_num(dataset) |> show_plot()

# Categorical columns
inspect_cat(dataset) |> show_plot()

table1(~ AF + age + PR | treatment, data = dataset)
