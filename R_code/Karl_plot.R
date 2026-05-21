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
library(showtext)
library(ggrepel)

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
    # age = as.integer(age),
    # age = ifelse(age >= 0 & age <= 120, age, NA),
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

snp_cols <- grep("^rs", names(dataset), value = TRUE)


############################# KARL TASK PLOT
extremes <- dataset %>%
  select(AF, sex, treatment, ID) %>%
  group_by(treatment, sex) %>%
  filter(AF == max(AF) | AF == min(AF)) %>%
  ungroup()

## Add font
font_add(
  family = "times",
  regular = here::here(
    "figs",
    "Times New Roman.ttf"
  )
)
# Set running
showtext_auto()

col_pal <- c("#006699", "#B35900")

# AF vs Treatment (by sex)
karl_plot_01 <- ggplot(dataset, aes(x = treatment, y = AF, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(
    data = extremes,
    aes(group = sex),
    shape = 21,
    size = 1,
    color = "black",
    position = position_dodge(width = 1)
  ) +
  geom_text(
    data = extremes,
    aes(
      label = paste0("ID: ", ID),
      group = sex
    ),
    position = position_dodge(width = 1),
    hjust = -0.1,
    size = 20,
    family = "times"
  ) +
  stat_summary(fun = "mean",
               geom = "point",
               shape = 22,
               size = 3,
               fill = "white",
               aes(group = sex),
               position = position_dodge(width = 1))+
  labs(title = "Boxplot of AF against treatment",
       subtitle = "White squares are sample means",
       x = "Treatment",
       y = "Activating factor",
       caption = "Produced by Jono and Joseph"
       )+
  theme(
    text = element_text(size = 70, family = "times"),

    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.line = element_line(colour = "black"),

    legend.background = element_blank(),
    legend.key = element_blank(),

    legend.position = "top",
    plot.caption = element_text(hjust = 1)
  ) +
  scale_fill_manual(values = col_pal)

ggsave(
  filename = here::here('figs', "Karl_plot_01.tiff"),
  plot = karl_plot_01,
  width = 9,
  height = 6,
  units = "in",
  dpi = 500,
  compression = "lzw"
)

karl_plot_01
