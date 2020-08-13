############################################################################## #
## ANALYSES ####################################################################
############################################################################## #


## PACKAGES AND FUNCTIONS ######################################################

library(tidyverse)
library(mgcv)
library(MuMIn)
library(stargazer)

# Set working directory
setwd("/home/danielschoenig/projects/forestchange_col_pa/")


## PARAMETERS ##################################################################

# Set forest cover threshold for analysis. This will assure that the correct
# files are loaded. For the reanalysis we used 1% and 50% cover thresholds.
fc_threshold <- 1


############################################################################## #
## PREPARE FOREST LOSS DATA FOR ANALYSIS #######################################
############################################################################## #

## LOAD DATASETS ###############################################################

forestchange_pa <- read_csv(paste0("results/forestchange/forestchange_pa_",
                                   str_pad(fc_threshold, 2, "left", "0"),
                                   ".csv"),
                            col_types = "icdddd")

forestchange_bz <- read_csv(paste0("results/forestchange/forestchange_bz_",
                                   str_pad(fc_threshold, 2, "left", "0"),
                                   ".csv"),
                            col_types = "icdddd")

forestchange_col <- read_csv(paste0("results/forestchange/forestchange_col_",
                                    str_pad(fc_threshold, 2, "left", "0"),
                                    ".csv"),
                             col_types = "dddd")

## CORRECT NAMES AND IDS #######################################################

# Remove non-terrestrial PAs and arrange by name
forestchange_pa <- forestchange_pa %>%
  filter(!pa_id %in% c(12, 16, 23, 41)) %>%
  arrange(nombre)
forestchange_bz <- forestchange_bz %>%
  filter(!pa_id %in% c(12, 16, 23, 41)) %>%
  arrange(nombre)

# Combine data on buffer zones and protected areas
forestchange <- bind_rows(pa = forestchange_pa,
                          bz = forestchange_bz,
                          .id = "type") %>%
         # Correct names
  mutate(name = str_replace_all(nombre,
                                c("Cahuinari" = "Cahuinarí",
                                  "Catatumbo Bari" = "Catatumbo Barí",
                                  "Volcanico Dona" = "Volcánico Doña",
                                  "Guacharos" = "Guácharos",
                                  "Farallones" = "Los Farallones",
                                  "Katios" = "Katíos",
                                  "Orquideas" = "Orquídeas",
                                  "Purace" = "Puracé",
                                  "Rio Pure" = "Río Puré",
                                  "Serrania" = "Serranía",
                                  "Tama" = "Tamá",
                                  "Tatama" = "Tatamá",
                                  "Yaigoje" = "Yaigojé",
                                  "Utria" = "Utría")),
         type = factor(type, levels = c("pa", "bz"))) %>%
  arrange(type, name) %>%
  mutate(pa_id = rep(1:(nrow(.)/2), 2)) %>%
  select(type, pa_id, name, area_forest_km2, no_loss_km2, 
         loss_1315_km2, loss_1618_km2)


# Generate labels for facets used in plotting
labs_type <- c(pa = "Protected areas", bz = "Buffer zones")



## CALCULATION OF CHANGE TRAJECTORIES ##########################################

# Relative change as in Clerici et al. 2020 for replication
forestchange_c1 <- forestchange %>%
  mutate(change_loss_km2 = loss_1618_km2 - loss_1315_km2,
         change_loss_percent = 100 * change_loss_km2 / loss_1315_km2)

# Per-period forest loss proportional to forest area, for extended analysis
forestchange_prop <- forestchange %>%
  mutate_at(vars(contains("_km2")), as.numeric) %>%
  pivot_longer(cols = starts_with("loss_"),
               names_to = "lossperiod",
               names_pattern = "loss_(.*)_km2",
               values_to = "loss_km2") %>%
  # Deforestation (% of forested area)
  mutate(loss_prop = 100 * loss_km2 / area_forest_km2)



## CALCULATION OF REFERENCE TRENDS #############################################

# Harmonize column order (same as `forestchange` data frame)
forestchange_col <- forestchange_col %>%
  select(!!names(select(forestchange, -type, -pa_id, -name)))

# Calculate forest loss totals per type and period
forestchange_sum <- forestchange %>%
  group_by(type) %>%
  select(-pa_id, -name) %>%
  summarise_all(sum) %>%
  arrange(type) %>%
  column_to_rownames("type") %>%
  as.matrix()

# Calculate reference values at the national level, outside PAs, and outside PAs
# with buffer areas
forestchange_ref <-
  do.call(rbind,
          list(national = as.matrix(forestchange_col),
               outside_pa = as.matrix(forestchange_col) -
                 forestchange_sum[2,],
               outside_pa_bz = as.matrix(forestchange_col) -
                 colSums(forestchange_sum))) %>%
  as_tibble() %>%
  mutate(scale = as.factor(c("national", "outside_pa", "outside_pa_bz"))) %>%
  select(scale, everything())
# Check that order of rows is correct (output should be `TRUE`):
is.unsorted(rev(forestchange_ref$area_forest_km2)) == FALSE

# Relative change at national level
forestchange_ref_c1 <- forestchange_ref %>%
  mutate(change_loss_percent =
           100 * (loss_1618_km2 - loss_1315_km2) / loss_1315_km2)

# Forest loss proportional to area at the national level
forestchange_ref_prop <- forestchange_ref %>%
  mutate_at(vars(contains("_km2")), as.numeric) %>%
  pivot_longer(cols = starts_with("loss_"),
               names_to = "lossperiod",
               names_pattern = "loss_(.*)_km2",
               values_to = "loss_km2") %>%
  # Deforestation (% of forested area)
  mutate(loss_prop = 100 * loss_km2 / area_forest_km2)



############################################################################## #
## COMPARISON OF FOREST LOSS DATA TO ORIGINAL ANALYSIS #########################
############################################################################## #

## COMPARISON OF FOREST LOSS PER PERIOD ########################################

# Load results from original study
forestchange_clerici <-
  read_csv("data/forestchange_Clerici_et_al_2020.csv", 
           col_types = "ficdddd") %>%
  arrange(type, pa_id)

data.frame(n1=forestchange_c1$name, n2=forestchange_clerici$name)

# Check for differences in percent of value from original analysis
forestchange.diff <- forestchange_c1 %>%
  select(-name) %>%
  left_join(forestchange_clerici,
             by = c("pa_id", "type"),
             suffix = c(".r", ".o")) %>%
  # Calculate absolute and relative (to original study) differences
  mutate(loss_1315_km2.diff = (loss_1315_km2.r - loss_1315_km2.o),
         loss_1618_km2.diff = (loss_1618_km2.r - loss_1618_km2.o),
         change_loss_km2.diff = change_loss_km2.r - change_loss_km2.o,
         change_loss_percent.diff = 
           change_loss_percent.r - change_loss_percent.o,
         loss_1315_km2.diff_rel = loss_1315_km2.diff / loss_1315_km2.o,
         loss_1618_km2.diff_rel = loss_1618_km2.diff / loss_1618_km2.o,
         change_loss_km2.diff_rel = change_loss_km2.diff / change_loss_km2.o,
         change_loss_percent.diff_rel =
           change_loss_percent.diff / change_loss_percent.o,
         ) %>%
  select(type, pa_id, name, area_forest_km2, no_loss_km2,
         loss_1315_km2.r, loss_1315_km2.o,
         loss_1315_km2.diff, loss_1315_km2.diff_rel,
         loss_1618_km2.r, loss_1618_km2.o,
         loss_1618_km2.diff, loss_1618_km2.diff_rel,
         change_loss_km2.r, change_loss_km2.o,
         change_loss_km2.diff, change_loss_km2.diff_rel,
         change_loss_percent.r, change_loss_percent.o,
         change_loss_percent.diff, change_loss_percent.diff_rel)

# Differences expressed as percent of value in original analysis
forestchange.diff_p <- forestchange.diff %>%
  mutate(loss_1315_km2.diff = loss_1315_km2.diff / loss_1315_km2.o,
         loss_1618_km2.diff = loss_1618_km2.diff / loss_1618_km2.o,
         change_loss_km2.diff = change_loss_km2.diff / change_loss_km2.o,
         change_loss_percent.diff =
           change_loss_percent.diff / change_loss_percent.o) %>%
  arrange(type, pa_id)


# Visual assessment
forestchange.diff %>%
  ggplot(aes(x = loss_1315_km2.r, y = loss_1315_km2.o)) +
  geom_point()

forestchange.diff %>%
  ggplot(aes(x = loss_1618_km2.r, y = loss_1618_km2.o)) +
  geom_point()

# Pearson's r for raw results
with(forestchange.diff, cor(loss_1315_km2.r, loss_1315_km2.o))
with(forestchange.diff, cor(loss_1618_km2.r, loss_1618_km2.o))

# Compare totals
apply(forestchange.diff[, c("loss_1315_km2.r", "loss_1315_km2.o",
                            "loss_1618_km2.r", "loss_1618_km2.o")],
      2, sum)

# Compare relative change (median and mean)
apply(forestchange.diff[,c("change_loss_percent.r",
                           "change_loss_percent.o")], 2, median)
apply(forestchange.diff[,c("change_loss_percent.r",
                           "change_loss_percent.o")], 2, mean)

## RELATIVE CHANGE IN DEFORESTATION ############################################

# Distribution of values
hist(forestchange_c1$change_loss_percent)

# Distribution is heavy-tailed and not symmetrical. Therefore, instead of the
# Wilcox test a sign test should be used.

change_loss_pa <- filter(forestchange_c1, type == "pa") %>%
  pull(change_loss_percent)

change_loss_bz <- filter(forestchange_c1, type == "bz") %>%
  pull(change_loss_percent)

# Is the relative change greater than 0? (1-sided)
binom.test(sum(change_loss_pa > 0),
           sum(change_loss_pa != 0),
           p = 0.5,
           alternative = "greater")

# Different to reference trend (outside pa and bz)? (2-sided)
binom.test(
  sum(change_loss_pa > with(forestchange_ref_c1,
                            change_loss_percent[scale == "outside_pa_bz"])),
  sum(change_loss_pa != with(forestchange_ref_c1,
                             change_loss_percent[scale == "outside_pa_bz"])),
  p = 0.5)

# Even though the Wilcox test is discouraged, it delivers qualitatively the same
# results.
# Is increase greater than 0?
filter(forestchange_c1, type == "pa") %>%
  pull(change_loss_percent) %>%
  wilcox.test(alternative = "greater")
# Different from reference trend
filter(forestchange_c1, type == "pa") %>%
  pull(change_loss_percent) %>%
  wilcox.test(mu = with(forestchange_ref_c1,
                        change_loss_percent[scale == "outside_pa_bz"]))

# Reproduce figure 1 of original study, adding reference trends
forestchange_c1 %>%
  mutate(change_loss_percent = as.numeric(change_loss_percent)) %>%
  ggplot(aes(y = change_loss_percent)) +
  facet_wrap("type", labeller = as_labeller(labs_type)) +
  ## geom_violin()
  geom_boxplot() +
  geom_hline(data = filter(forestchange_ref_c1,
                           scale != "outside_pa"),
             aes(yintercept = change_loss_percent, col = scale),
             size = 0.5) +
  scale_y_continuous(breaks = seq(-200, 600, 200),
                     # minor_breaks = NULL,
                     limits = c(NA, 700)) +
  scale_x_discrete(breaks = c("pa", "bz"),
                   labels = c("Protected areas", "Buffer zones")) +
  scale_color_discrete(breaks = c("national", "outside_pa_bz"),
                       labels = c("National",
                                  str_wrap("Outside protected areas
                                           and buffer zones", 20))) +
  labs(x = NULL,
       y = "Change in deforested area (percent)",
       col = "Reference trend") +
  theme_bw()


############################################################################## #
## EXTENDED ANALYSIS ###########################################################
############################################################################## #

# In the following, forest loss is expressed as a proportion of total forested
# area. This allows to fully leverage the data and it is more adequate than
# using relative changes only: the same relative change leads to very different
# trajectories depending on the initial deforestation rate. This allows forest
# loss to be compared both between periods and against reference trends.


## INITIAL VISUAL ASSESSMENT ###################################################

# First assessment in the style of Fig. 1 of the original study, but using
# proportional forest loss.
forestchange_prop %>%
  ggplot(aes(y = loss_prop, x = lossperiod)) +
  facet_wrap("type", labeller = as_labeller(labs_type)) +
  geom_boxplot(alpha = 0.8) +
  geom_line(data =
              ## forestchange_ref_prop,
              filter(forestchange_ref_prop,
                     scale != "outside_pa"),
            aes(x = lossperiod, y = loss_prop,
                group = scale, col = scale),
            size = 1.2) +
  scale_color_discrete(breaks = c("national", "outside_pa_bz"),
                       labels = c("National",
                                  str_wrap("Outside protected areas
                                           and buffer zones", 20))) +
  scale_x_discrete(breaks = c("1315", "1618"),
                   labels = c("before", "after")) +
  ylim(NA, 3) +
  labs(x = "Peace agreement",
       y = "% of forest area lost",
       col = "Reference trend") +
  theme_bw()

# Identify individual PAs
forestchange_prop %>%
  ggplot(aes(y = loss_prop, x = lossperiod)) +
  facet_wrap("type", labeller = as_labeller(labs_type)) +
  geom_line(aes(group = pa_id, col = pa_id),
            alpha = 0.3) +
  geom_text(aes(label = pa_id)) +
  ## geom_boxplot(alpha = 0.8) +
  scale_x_discrete(breaks = c("1315", "1618"),
                   labels = c("before", "after")) +
  labs(title = "Forest loss",
       x = "Peace agreement",
       y = "% of forst area lost",
       col = "Protected Area")

  # PA 37 is a strong outlier


## REGRESSION MODEL ############################################################

# Heavy-tailed and asymmetrical distribution with concrete probability mass at 0
# (see histogram above) can be accomodated best by a Tweedie distribution.

# The following models are fitted with the GAM function of mgcv, although all
# predictors are strictly parametric (i.e. no smooth terms are used) and only
# random effects are fitted as penalized terms.

# Intercept only
mg1 <- gam(loss_prop ~ 1,
           family = tw,
           method = "REML",
           data = forestchange_prop)
summary(mg1)

# Main effects without random effect
mg2 <- gam(loss_prop ~ lossperiod + type,
           family = tw,
           method = "REML",
           data = forestchange_prop)
summary(mg2)

# Main effects with random effect
mg3 <- gam(loss_prop ~ lossperiod + type + s(pa_id, bs = "re"),
           family = tw,
           method = "REML",
           data = forestchange_prop)
summary(mg3)

# Main effects with interaction and random effect
mg4 <- gam(loss_prop ~ lossperiod * type + s(pa_id, bs = "re"),
           family = tw,
           method = "REML",
           data = forestchange_prop)
summary(mg4)

AICc(mg1, mg2, mg3, mg4)
# Third model appears best

# Model checking plots and diagnostics
gam.check(mg3)
# Nothing of great concern

# Test significance of terms
anova(mg3)

# Variance componentes
gam.vcomp(mg3)

# Removing the strongest outlier does not affect the model much
mg3_37rm <- gam(loss_prop ~ lossperiod + type + s(pa_id, bs = "re", k = 50),
                family = tw(link = log),
                data = filter(forestchange_prop, pa_id != 37))
summary(mg3_37rm)

## Model results and comparison to reference trend and counterfactual ##########

# Model predictions on response scale
forestchange_prop_pred <-
  expand.grid(lossperiod = unique(forestchange_prop$lossperiod),
              type = unique(forestchange_prop$type),
              pa_id = 0)
predicted_loss_prop <-
  predict(mg3, expand.grid(lossperiod = unique(forestchange_prop$lossperiod),
                      type = unique(forestchange_prop$type),
                      pa_id = 0),
          type = "link",
          se.fit = TRUE)
forestchange_prop_pred <- forestchange_prop_pred %>%
  mutate(loss_prop = exp(predicted_loss_prop$fit),
         loss_prop_ci_l =
           exp(predicted_loss_prop$fit - 1.96*predicted_loss_prop$se.fit),
         loss_prop_ci_u =
           exp(predicted_loss_prop$fit + 1.96*predicted_loss_prop$se.fit)) %>%
  select(-pa_id)

# Prepare reference trend and counterfactual
counterfactual <- forestchange_prop_pred %>%
  select(lossperiod, type, loss_prop) %>%
  arrange(lossperiod, type)
loss_prop_outside_pa_bz <- forestchange_ref_prop %>%
  filter(scale == "outside_pa_bz") %>%
  select(lossperiod, loss_prop) %>%
  arrange(lossperiod)
counterfactual$loss_prop[c(3,4)] <- counterfactual$loss_prop[c(1,2)] +
    loss_prop_outside_pa_bz$loss_prop[2] - loss_prop_outside_pa_bz$loss_prop[1]

# Visual comparison
plot_forestchange_model <- forestchange_prop_pred %>%
  ggplot(aes(y = loss_prop, x = lossperiod, group = type)) +
  facet_wrap("type", labeller = as_labeller(labs_type)) +
  geom_point(data = forestchange_prop,
             col = "grey20", size = 0.5) +
  geom_line(data = filter(forestchange_ref_prop,
                          !scale %in% c("outside_pa", "national")),
            aes(x = lossperiod, y = loss_prop,
                group = scale),
            size = 0.5, col = "#DF536B",
            linetype = "dashed"
            ) +
  geom_line(data = counterfactual,
            aes(x = lossperiod, y = loss_prop,
                group = type),
            size = 0.5, col = "#2297E6",
            linetype = "dashed"
            ) +
  geom_ribbon(aes(ymin = loss_prop_ci_l, ymax = loss_prop_ci_u),
              alpha = 0.15) +
  geom_line(size = 0.7) +
  scale_x_discrete(breaks = c("1315", "1618"),
                   labels = c("before", "after")) +
  scale_color_discrete(breaks = c("national", "outside_pa_bz"),
                       labels = c("national",
                                  str_wrap("outside protected areas
                                           and buffer zones", 20))) +
  ylim(NA, 4.5) +
  labs(x = "Peace agreement",
       y = "% of forest area lost",
       col = "Reference trend") +
  theme_classic(base_family = "IBMPlexSans") +
  theme(strip.text.x = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 11, face = "bold"),
        strip.background =
          element_rect(color=NA, fill="grey90", size=1.5, linetype="solid"))

plot_forestchange_model

# Not plotted:
filter(forestchange_prop, loss_prop > 4.5) %>% arrange(type, lossperiod)
paste0("results/plots/raw/forestchange_model_",
       str_pad(fc_threshold, 2, "left", "0"),
       ".svg")
ggsave(file = paste0("results/plots/raw/forestchange_model_",
                     str_pad(fc_threshold, 2, "left", "0"),
                     ".svg"),
       plot = plot_forestchange_model,
       width = 15, height = 12, units = "cm", dpi = 300)

## EXTRACT QUANTITATIVE INFORMATION FOR MANUSCRIPT ##############################

# Comparison of model results
comparisons <- forestchange_prop_pred %>%
  left_join(loss_prop_outside_pa_bz,
            by = "lossperiod",
            suffix = c("", "_outside_pa_bz")) %>%
  mutate(vs_out = loss_prop - loss_prop_outside_pa_bz,
         vs_out_ci_l = loss_prop_ci_l - loss_prop_outside_pa_bz,
         vs_out_ci_u = loss_prop_ci_u - loss_prop_outside_pa_bz) %>%
  left_join(counterfactual,
            by = c("lossperiod", "type"),
            suffix = c("", "_counterfactual")) %>%
  mutate(vs_cf = loss_prop - loss_prop_counterfactual,
         vs_cf_ci_l = loss_prop_ci_l - loss_prop_counterfactual,
         vs_cf_ci_u = loss_prop_ci_u - loss_prop_counterfactual) %>%
  select(-loss_prop_outside_pa_bz, -loss_prop_counterfactual) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  mutate(loss_prop = str_c(loss_prop, " (", loss_prop_ci_l, ", ", loss_prop_ci_u, ")"),
         vs_out = str_c(vs_out, " (", vs_out_ci_l, ", ", vs_out_ci_u, ")"),
         vs_cf =  str_c(vs_cf, " (", vs_cf_ci_l, ", ", vs_cf_ci_u, ")")) %>%
  select(type, lossperiod, loss_prop, vs_out, vs_cf)

# Model effects
mg3_vcomp <- gam.vcomp(mg3)
mg3_sum <- summary(mg3)
effects <- mg3_sum$p.table[,1:2] %>%
  as.data.frame() %>%
  rownames_to_column(var = "parameter") %>%
  rename(f.coef = Estimate,
         f.se = "Std. Error") %>%
  mutate(f.ci_l = f.coef - 1.96 * f.se,
         f.ci_u = f.coef + 1.96 * f.se,
         r.sd = NA,
         r.ci_l = NA,
         r.ci_u = NA) %>%
  add_case(parameter = row.names(mg3_vcomp)[1],
           r.sd = mg3_vcomp[1,1],
           r.ci_l = mg3_vcomp[1,2],
           r.ci_u = mg3_vcomp[1,3]) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  mutate(f.ci = str_c(f.ci_l, ", ", f.ci_u),
         r.ci = str_c(r.ci_l, ", ", r.ci_u)) %>%
  select(parameter, f.coef, f.se, f.ci, r.sd, r.ci)

# Overview of loss figures
forestchange %>%
  select(-no_loss_km2) %>%
  pivot_wider(names_from = type, 
              values_from = c(area_forest_km2, loss_1315_km2, loss_1618_km2),
              names_sep = ".") %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  select(name, 
         pa_id, 
         area_forest_km2.pa, 
         loss_1315_km2.pa,
         loss_1618_km2.pa,
         area_forest_km2.bz,
         loss_1315_km2.bz,
         loss_1618_km2.bz) %>%
  stargazer(summary = FALSE, rownames = FALSE)
