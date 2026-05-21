pacman::p_load(
  readxl,
  tidymodels,
  tidyverse,
  gt
)

# Loading
dataset1 <-
  read_excel(here::here("raw-data", "FITFR-final-data-v3.xlsx"),
             skip = 1,
             n_max = 100)

dataset2 <-
  read_excel(here::here("raw-data", "FITFR-extra-data.xlsx"),
             n_max = 10)

dataset <- rbind(dataset1, dataset2)

# Recoding
dataset <- dataset %>%
  mutate(
    age = as.integer(age),
    age = ifelse(age >= 0 & age <= 120, age, NA),
    sex = ifelse(tolower(sex) %in% c("male", "m"),
                 "Male",
                 ifelse(tolower(sex) %in% c("female", "f"), "Female", NA)),
    treatment = toupper(treatment)
    )

# Remove rows where there are missing values
dataset <- dataset %>%
  filter(complete.cases(.))

# Convert categorical columns to factors
snp_cols <- grep("^rs", names(dataset), value = TRUE)
dataset[snp_cols] <- lapply(dataset[snp_cols], as.factor)
dataset$sex <- as.factor(dataset$sex)
dataset$treatment <- as.factor(dataset$treatment)
dataset$smoker <- as.factor(dataset$smoker)
dataset$PR <- as.factor(dataset$PR)

# Remove columns where all entries are the same
to_remove <- dataset %>%
  select(where(~ length(unique(.)) < 2))
dataset <- dataset %>%
  select(where(~ length(unique(.)) >= 2)) %>%
  select(-ID)

# PCA
snp_data <- dataset[, snp_cols]
snp_numeric <- as.data.frame(lapply(snp_data, function(x) {
  as.numeric(as.factor(x)) - 1
}))
pca <- prcomp(snp_numeric, center = TRUE, scale. = TRUE)
pc_scores <- as.data.frame(pca$x[, 1:2])  # first 5 PCs
colnames(pc_scores) <- paste0("PC", 1:2)
df <- dataset %>%
  select(-snp_cols) %>%
  cbind(pc_scores)
model <- glm(PR ~ age + sex, data = df, family = binomial)
summary(model)

# Setting up cross validation
num_folds <- 5
dataset_cv <- vfold_cv(dataset, v = num_folds)


######################################## AF No Interaction
AF_recipeNo <-
  recipe(AF ~ ., data = dataset) |>
  step_rm(PR) |>
  step_dummy(all_nominal_predictors()) |>
  step_normalize(all_numeric_predictors())

AF_modelNo <- linear_reg(penalty = tune(), mixture = 1) |>
  set_mode("regression") |>
  set_engine("glmnet")

AF_wfNo <- workflow(AF_recipeNo, AF_modelNo)

doParallel::registerDoParallel()

search_grid <- grid_regular(penalty(), levels = 50)

AF_tuneNo <- tune_grid(
  AF_wfNo,
  resamples = dataset_cv,
  grid = search_grid
)

collect_metrics(AF_tuneNo)

AF_penaltyNo <- select_best(AF_tuneNo, metric = "rmse")

AF_wfNo <- AF_wfNo |>
  finalize_workflow(AF_penaltyNo)

AF_fitNo <- AF_wfNo |> fit(dataset)
AF_fitNo |>
  extract_fit_parsnip() |>
  vip::vi() |>
  filter(Importance > 0.5) |>
  mutate(
    Variable = fct_reorder(Variable, Importance)
  ) |>
  ggplot(aes(Importance, Variable, fill = Sign)) +
  geom_col() +
  harrypotter::scale_fill_hp_d("Ravenclaw")

af_augNo <- dataset |>
  mutate(
    .fitted = predict(AF_fitNo, new_data = dataset)$.pred,
    .resid = AF - .fitted
  )

ggplot(af_augNo, aes(.fitted, .resid, colour = treatment)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  labs(
    x = "Fitted values",
    y = "Residuals",
    colour = "Treatment",
    title = "AF model: Residuals vs Fitted"
  )

####################################### AF
AF_recipe <-
  recipe(AF ~ ., data = dataset) |>
  step_rm(PR) |>
  step_dummy(all_nominal_predictors()) |>
  step_interact(terms = ~ starts_with("treatment"):age) |>
  step_normalize(all_numeric_predictors())

AF_model <- linear_reg(penalty = tune(), mixture = 1) |>
  set_mode("regression") |>
  set_engine("glmnet")

AF_wf <- workflow(AF_recipe, AF_model)

doParallel::registerDoParallel()

search_grid <- grid_regular(penalty(), levels = 50)

AF_tune <- tune_grid(
  AF_wf,
  resamples = dataset_cv,
  grid = search_grid
)

collect_metrics(AF_tune)
AF_tune |> autoplot()

AF_summary <- show_best(AF_tune, metric = "rmse") |>
  mutate(
    sd = std_err * sqrt(num_folds),
    relative_sd = sd / mean
  ) |>
  transmute(
    penalty,
    metric = .metric,
    mean = round(3),
    sd = round(sd, 3),
  )
AF_summary |>
  gt() |>
  fmt_scientific(
    columns = penalty,
    decimals = 2
  )

AF_penalty <- select_best(AF_tune, metric = "rmse")


AF_wf <- AF_wf |>
  finalize_workflow(AF_penalty)


AF_fit <- AF_wf |> fit(dataset)
AF_fit |>
  extract_fit_parsnip() |>
  vip::vi() |>
  filter(Importance > 0.001) |>
  mutate(
    Variable = fct_reorder(Variable, Importance)
  ) |>
  ggplot(aes(Importance, Variable, fill = Sign)) +
  geom_col() +
  harrypotter::scale_fill_hp_d("Ravenclaw")

af_aug <- dataset |>
  mutate(
    .fitted = predict(AF_fit, new_data = dataset)$.pred,
    .resid = AF - .fitted
  )

ggplot(af_aug, aes(.fitted, .resid, colour = treatment)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  labs(
    x = "Fitted values",
    y = "Residuals",
    colour = "Treatment",
    title = "AF model: Residuals vs Fitted"
  )

coef_table_AF <- AF_fit |>
  extract_fit_parsnip() |>
  tidy()

nonzero_coef_AF <- coef_table_AF |>
  filter(estimate != 0)

##################################################### PR
PR_recipe <-
  recipe(PR ~ ., data = dataset) |>
  step_rm(AF) |>
  step_dummy(all_nominal_predictors())|>
  step_zv(all_numeric_predictors()) |>
  step_normalize(all_numeric_predictors())

PR_model <- logistic_reg(penalty = tune(), mixture = tune()) |>
  set_mode("classification") |>
  set_engine("glmnet")

PR_wf <- workflow(PR_recipe, PR_model)

doParallel::registerDoParallel()

search_grid <- grid_space_filling(
  penalty(),
  mixture(),
  size = 50
)

PR_tune <- tune_grid(
  PR_wf,
  resamples = dataset_cv,
  grid = search_grid
)

PR_tune |> autoplot()

PR_summary <- show_best(PR_tune, metric = "brier_class") |>
  mutate(
    sd = std_err * sqrt(num_folds),
    relative_sd = sd / mean
  ) |>
  transmute(
    penalty,
    metric = .metric,
    mean = round(mean, 3),
    sd = round(sd, 3),
    relative_sd = round(relative_sd, 3)
  )
PR_summary

PR_penalty <- select_best(PR_tune, metric = "brier_class")

PR_wf <- PR_wf |>
  finalize_workflow(PR_penalty)


PR_fit <- PR_wf |> fit(dataset)
PR_fit |>
  extract_fit_parsnip() |>
  vip::vi() |>
  filter(Importance > 1) |>
  mutate(
    Variable = fct_reorder(Variable, Importance)
  ) |>
  ggplot(aes(Importance, Variable, fill = Sign)) +
  geom_col() +
  harrypotter::scale_fill_hp_d("Ravenclaw")


# Calibration curve
pr_probs <- predict(PR_fit, dataset, type = "prob") |>
  bind_cols(dataset |> select(PR))

cal_data <- pr_probs |>
  mutate(bin = ntile(.pred_yes, 10)) |>
  group_by(bin) |>
  summarise(
    mean_pred = mean(.pred_yes),
    observed = mean(as.numeric(PR == "yes"))
  )

ggplot(cal_data, aes(mean_pred, observed)) +
  geom_point() +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Predicted probability",
    y = "Observed proportion",
    title = "Calibration curve (PR model)"
  )

coef_table_PR <- PR_fit |>
  extract_fit_parsnip() |>
  tidy()

nonzero_coef_PR <- coef_table_PR |>
  filter(estimate != 0)
nonzero_coef_PR |>
  gt() |>
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  )
