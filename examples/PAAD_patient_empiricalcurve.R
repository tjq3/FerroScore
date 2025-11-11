######Ferroptosis clinical relevance in PAAD patients
######################################################



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("tidyverse", "readr","ggplot2","mgcv","rstatix",
                       "PMCMRplus","vcd","ggpubr","coin","reshape2",
                       "broom","dplyr","glmnet","caret","ggrepel",
                       "patchwork","clusterProfiler","org.Hs.eg.db",
                       "GSVA","GSEABase","msigdbr","lightgbm","minpack.lm"))

library(tidyverse)
library(readr)
library(ggplot2)
library(mgcv)
library(rstatix)
library(PMCMRplus)
library(vcd)       
library(ggpubr)
library(coin)
library(reshape2)
library(broom)
library(dplyr)
library(glmnet)       
library(caret)        
library(ggrepel)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)         
library(GSEABase)   
library(msigdbr)  
library(lightgbm)     
library(minpack.lm)


output_dir <- "my_results/empirical_curve"


#Ferroptosis Score & Survival Time empirical curve
ferro_data <- read_csv("my_data/ferroptosis_scores_PAAD.csv")
ferro_scores <- ferro_data[4,]
ferro_scores_transposed <- ferro_scores %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "patient_id") %>% 
  mutate(patient_id = gsub("\\.", "-", patient_id)) %>%
  mutate(patient_id = substr(patient_id, 1, 12))
ferro_scores<-ferro_scores_transposed[-1,]
colnames(ferro_scores)[2] <- "Ferroptosis_score"

clinical_data <- read_csv("my_data/TCGA_PAAD_clinical.csv")
survival_time <- clinical_data %>% 
  select(patient_id = 1, survival_time = 16)

merged_data_survival <- inner_join(survival_time, ferro_scores, by = "patient_id")
write.csv(merged_data_survival, file = file.path(output_dir, "merged_data_PAAD_survival.csv"), row.names = TRUE)




#Using LOESS fitting
merged_data_survival$Ferroptosis_score <- as.numeric(merged_data_survival$Ferroptosis_score)
merged_data_survival$survival_time <- as.numeric(merged_data_survival$survival_time)
data_clean <- na.omit(merged_data_survival[, c("Ferroptosis_score", "survival_time")])
span_range <- seq(0.1, 1.0, by = 0.05)
gcv_values <- numeric(length(span_range))
trace_hat_values <- numeric(length(span_range))

#Loop to calculate GCV for each span
for (i in seq_along(span_range)) {
  span_val <- span_range[i]
  loess_fit <- loess(survival_time ~ Ferroptosis_score, 
                     data = data_clean,
                     span = span_val,
                     degree = 2,
                     family = "gaussian")
  n <- nrow(data_clean)
  sse <- sum(loess_fit$residuals^2)
  trace_hat <- loess_fit$trace.hat
  gcv <- (sse / n) / ((1 - trace_hat / n)^2)
  gcv_values[i] <- gcv
  trace_hat_values[i] <- trace_hat
  cat(sprintf("%.2f\t%.4f\t\t%.2f\n", span_val, gcv, trace_hat))
}

#Find optimal span corresponding to minimum GCV
optimal_idx <- which.min(gcv_values)
optimal_span <- span_range[optimal_idx]
min_gcv <- gcv_values[optimal_idx]

par(mar = c(5, 5, 4, 2) + 0.1)
plot(span_range, gcv_values, type = "b", pch = 16, col = "blue", lwd = 2,
     xlab = "Span Value", ylab = "GCV Value", cex.lab = 1.2, cex.axis = 1.1,
     main = "Optimal Span Selection Based on GCV Criterion", cex.main = 1.3)
abline(v = optimal_span, col = "red", lty = 2, lwd = 2)
points(optimal_span, min_gcv, col = "red", pch = 16, cex = 2)
text(optimal_span, min_gcv, 
     labels = sprintf("Optimal span: %.2f", optimal_span),
     pos = 3, col = "red", cex = 1.1)
grid()

#Refit final model using optimal span
pdf(file.path(output_dir, "FerroScores&SurvivalTime_empiricalcurve.pdf"), width = 9, height = 6)
ggplot(data_clean, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  geom_smooth(data = data_clean,
              method = "loess", 
              se = TRUE, 
              color = "#7291c0", 
              fill = "#b8c8df", 
              alpha = 0.3,
              linewidth = 1.2,
              span = optimal_span) +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.margin = margin(15, 15, 15, 15)
  )
dev.off()



#Fit y=f(u) analytical expression
#Find key segmentation points
ggplot(data_clean, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  geom_smooth(data = data_clean,
              method = "loess", 
              se = TRUE, 
              color = "#7291c0", 
              fill = "#b8c8df", 
              alpha = 0.3,
              linewidth = 1.2,
              span = optimal_span) +
  scale_x_continuous(breaks = seq(floor(min(data_clean$Ferroptosis_score)), 
                                  ceiling(max(data_clean$Ferroptosis_score)), 
                                  by = 0.1),
                     minor_breaks = seq(floor(min(data_clean$Ferroptosis_score)), 
                                        ceiling(max(data_clean$Ferroptosis_score)), 
                                        by = 0.01)) + 
  scale_y_continuous(breaks = seq(0, 
                                  ceiling(max(data_clean$survival_time)), 
                                  by = 500),
                     minor_breaks = seq(0, 
                                        ceiling(max(data_clean$survival_time)), 
                                        by = 100)) +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.3), 
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.margin = margin(15, 15, 15, 15)
  )

#Middle quadratic U-shaped part
data_u_shape <- data_clean[data_clean$Ferroptosis_score > 0.45 & data_clean$Ferroptosis_score <= 0.77, ]
model_quadratic <- nlsLM(survival_time ~ a*Ferroptosis_score^2 + b*Ferroptosis_score + c,
                         data = data_u_shape,
                         start = list(a = -1000, b = 1500, c = 1000),
                         control = nls.lm.control(maxiter = 1000))
params <- coef(model_quadratic)
a <- params["a"]
b <- params["b"]
c <- params["c"]
cat("Approximate expression for U-shaped curve (0.45 < x ≤ 0.77):\n")
cat(sprintf("y = %.2f*x² + %.2f*x + %.2f\n", a, b, c))
x_values <- seq(min(data_u_shape$Ferroptosis_score), 
                max(data_u_shape$Ferroptosis_score), 
                length.out = 100)
y_pred <- a * x_values^2 + b * x_values + c
u_curve_data <- data.frame(Ferroptosis_score = x_values, 
                           survival_time = y_pred)
pdf(file.path(output_dir, "y=f(u)_approximation_secondsegment.pdf"), width = 9, height = 6)
ggplot(data_clean, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  geom_smooth(data = data_clean,
              method = "loess", 
              se = TRUE, 
              color = "#7291c0", 
              fill = "#b8c8df", 
              alpha = 0.3,
              linewidth = 1.2,
              span = optimal_span) +
  geom_line(data = u_curve_data, 
            aes(x = Ferroptosis_score, y = survival_time),
            color = "red", 
            linewidth = 1.5, 
            linetype = "solid") +
  geom_vline(xintercept = 0.45, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  geom_vline(xintercept = 0.77, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.margin = margin(15, 15, 15, 15)
  )
dev.off()

#Front linear decreasing part
data_linear <- data_clean[data_clean$Ferroptosis_score <= 0.45, ]
model_linear <- lm(survival_time ~ Ferroptosis_score, data = data_linear)
params_linear <- coef(model_linear)
intercept <- params_linear["(Intercept)"]
slope <- params_linear["Ferroptosis_score"]
cat("Approximate expression for linear decreasing curve (x ≤ 0.45):\n")
cat(sprintf("y = %.2f + %.2f*x\n", intercept, slope))
x_values_linear <- seq(min(data_linear$Ferroptosis_score), 
                       max(data_linear$Ferroptosis_score), 
                       length.out = 100)
y_pred_linear <- intercept + slope * x_values_linear
linear_curve_data <- data.frame(Ferroptosis_score = x_values_linear, 
                                survival_time = y_pred_linear)
pdf(file.path(output_dir, "y=f(u)_approximation_firstsegment.pdf"), width = 9, height = 6)
ggplot(data_clean, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  geom_smooth(data = data_clean,
              method = "loess", 
              se = TRUE, 
              color = "#7291c0", 
              fill = "#b8c8df", 
              alpha = 0.3,
              linewidth = 1.2,
              span = optimal_span) +
  geom_line(data = u_curve_data, 
            aes(x = Ferroptosis_score, y = survival_time),
            color = "red", 
            linewidth = 1.5, 
            linetype = "solid") +
  geom_line(data = linear_curve_data, 
            aes(x = Ferroptosis_score, y = survival_time),
            color = "blue", 
            linewidth = 1.5, 
            linetype = "solid") +
  geom_vline(xintercept = 0.45, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  geom_vline(xintercept = 0.77, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.margin = margin(15, 15, 15, 15)
  )
dev.off()

#Back linear increasing part
data_linear_up <- data_clean[data_clean$Ferroptosis_score > 0.77, ]
model_linear_up <- lm(survival_time ~ Ferroptosis_score, data = data_linear_up)
params_linear_up <- coef(model_linear_up)
intercept_up <- params_linear_up["(Intercept)"]
slope_up <- params_linear_up["Ferroptosis_score"]
cat("Approximate expression for linear increasing curve (x > 0.77):\n")
cat(sprintf("y = %.2f + %.2f*x\n", intercept_up, slope_up))
x_values_up <- seq(min(data_linear_up$Ferroptosis_score), 
                   max(data_linear_up$Ferroptosis_score), 
                   length.out = 100)
y_pred_up <- intercept_up + slope_up * x_values_up
linear_up_data <- data.frame(Ferroptosis_score = x_values_up, 
                             survival_time = y_pred_up)
pdf(file.path(output_dir, "y=f(u)_approximation_thirdsegment.pdf"), width = 9, height = 6)
ggplot(data_clean, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  geom_smooth(data = data_clean,
              method = "loess", 
              se = TRUE, 
              color = "#7291c0", 
              fill = "#b8c8df", 
              alpha = 0.3,
              linewidth = 1.2,
              span = optimal_span) +
  geom_line(data = linear_curve_data, 
            aes(x = Ferroptosis_score, y = survival_time),
            color = "blue", 
            linewidth = 1.5, 
            linetype = "solid") +
  geom_line(data = u_curve_data, 
            aes(x = Ferroptosis_score, y = survival_time),
            color = "red", 
            linewidth = 1.5, 
            linetype = "solid") +
  geom_line(data = linear_up_data, 
            aes(x = Ferroptosis_score, y = survival_time),
            color = "green4", 
            linewidth = 1.5, 
            linetype = "solid") +
  geom_vline(xintercept = 0.45, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  geom_vline(xintercept = 0.77, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.margin = margin(15, 15, 15, 15)
  )
dev.off()




#Fit u=g(x) analytical expression
set.seed(123) 
#Load data
ferro_data <- read.csv("my_data/ferroptosis_scores_PAAD.csv", header = TRUE)
u_data <- data.frame(
  Sample = colnames(ferro_data)[-1],  
  u_score = as.numeric(ferro_data[4, -1])  
)
tpm_data <- read.csv("my_data/TCGA_PAAD_TPM.csv", header = TRUE, row.names = 1)

u_data$Sample <- gsub("\\.", "-", u_data$Sample)
colnames(tpm_data) <- gsub("\\.", "-", colnames(tpm_data))

common_samples <- intersect(u_data$Sample, colnames(tpm_data))
u_data_filtered <- u_data[u_data$Sample %in% common_samples, ]
tpm_data_filtered <- tpm_data[, common_samples]

u_data_filtered <- u_data_filtered[match(common_samples, u_data_filtered$Sample), ]
tpm_data_filtered <- tpm_data_filtered[, common_samples]

expr_matrix <- as.matrix(tpm_data_filtered)

#Prepare MSigDB pathway gene sets
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

h_gene_sets_clean <- h_gene_sets %>%
  group_by(gs_name, gene_symbol) %>%
  summarise(.groups = "drop")  

gene_set_list <- split(h_gene_sets_clean$gene_symbol, h_gene_sets_clean$gs_name)
gene_sets <- GeneSetCollection(lapply(names(gene_set_list), function(set_name) {
  unique_genes <- unique(gene_set_list[[set_name]])
  GeneSet(unique_genes, 
          setName = set_name, 
          geneIdType = SymbolIdentifier())
}))

#GSVA calculation for pathway activity scores
gsva_results <- gsva(expr_matrix, 
                     gene_sets,
                     method = "gsva",   
                     kcdf = "Gaussian", 
                     verbose = TRUE,
                     parallel.sz = 1) 
pathway_activity <- as.data.frame(t(gsva_results))
pathway_activity$Sample <- rownames(pathway_activity)

final_data <- merge(u_data_filtered, pathway_activity, by = "Sample")
rownames(final_data) <- final_data$Sample
final_data <- final_data[, -1]

#Data splitting
X <- final_data[, -1]  
y <- final_data$u_score  

train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, ]
X_test <- X[-train_index, ]
y_train <- y[train_index]
y_test <- y[-train_index]

#Elastic Net Regression
#Remove low variance pathways
low_var <- nearZeroVar(X_train, freqCut = 95/5)
if(length(low_var) > 0) {
  X_train_processed <- X_train[, -low_var]
  X_test_processed <- X_test[, -low_var]
  cat("Removed", length(low_var), "low variance pathways\n")
} else {
  X_train_processed <- X_train
  X_test_processed <- X_test
  cat("No low variance pathways to remove\n")
}
#Remove highly correlated pathways
cor_matrix <- cor(X_train_processed)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.9)
if(length(high_cor) > 0) {
  X_train_processed <- X_train_processed[, -high_cor]
  X_test_processed <- X_test_processed[, -high_cor]
  cat("Removed", length(high_cor), "highly correlated pathways\n")
} else {
  cat("No highly correlated pathways to remove\n")
}

set.seed(123)
enet_grid <- expand.grid(
  alpha = seq(0, 1, by = 0.1),  
  lambda = 10^seq(-4, 1, length.out = 50)  
)
enet_control <- trainControl(
  method = "cv", 
  number = 10,  
  allowParallel = FALSE,
  verboseIter = TRUE
)
enet_model_processed <- train(
  x = X_train_processed, 
  y = y_train,
  method = "glmnet",
  tuneGrid = enet_grid,
  trControl = enet_control,
  preProcess = c("center", "scale"),
  metric = "RMSE"
)
final_enet_model <- glmnet(as.matrix(X_train), y_train, 
                           alpha = 1, 
                           lambda = 0.0356,
                           standardize = TRUE) 
enet_coef <- as.matrix(coef(final_enet_model))
non_zero_indices <- which(enet_coef != 0)
selected_pathways <- rownames(enet_coef)[non_zero_indices]
selected_pathways <- selected_pathways[selected_pathways != "(Intercept)"]

cat("Number of pathways selected by Elastic Net:", length(selected_pathways), "\n")
if(length(selected_pathways) > 0) {
  cat("Selected pathways:", paste(selected_pathways, collapse = ", "), "\n")
}

#LightGBM
set.seed(123)
lgb_train <- lgb.Dataset(data = as.matrix(X_train), label = y_train)
lgb_test <- lgb.Dataset(data = as.matrix(X_test), label = y_test, free_raw_data = FALSE)

lgb_params <- list(
  objective = "regression",
  metric = "l2",
  boosting_type = "gbdt",
  learning_rate = 0.1,
  num_leaves = 31,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 5,
  lambda_l1 = 0.1,
  lambda_l2 = 0.1
)

lgb_model <- lgb.train(
  params = lgb_params,
  data = lgb_train,
  valids = list(test = lgb_test),
  nrounds = 100,
  early_stopping_rounds = 20,
  verbose = 1
)

lgb_importance <- lgb.importance(lgb_model, percentage = TRUE)
lgb_importance_df <- as.data.frame(lgb_importance)
lgb_importance_df$Pathway <- lgb_importance_df$Feature
lgb_importance_df <- lgb_importance_df[order(-lgb_importance_df$`Gain`), ]

#Visualize LightGBM feature importance (Top 15)
pdf(file.path(output_dir, "u=g(x)_pathway_LightGBM.pdf"), width = 9, height = 6)
p_lgb_importance <- ggplot(head(lgb_importance_df, 15), 
                           aes(x = reorder(Pathway, `Gain`), y = `Gain`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 15 Important Pathways (LightGBM)",
       x = "Pathway", y = "Importance (Gain)") +
  theme_minimal()
print(p_lgb_importance)
dev.off()




#Ferroptosis Score & Tumor Stage
clinical_data <- read_csv("my_data/TCGA_PAAD_clinical.csv")
stage_event <- clinical_data %>% 
  dplyr::select(patient_id = 1, stage_event = 11)

merged_data_stage_event <- inner_join(stage_event, ferro_scores, by = "patient_id")
write.csv(merged_data_stage_event, file = file.path(output_dir, "merged_data_PAAD_stage_event.csv"),row.names = TRUE)
merged_data_stagepre <- read.csv("my_data/merged_data_PAAD_stagepre.csv")

analysis_data_ss <- merged_data_stagepre %>%
  mutate(
    stage_event = factor(stage_event, 
                         levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage X"),
                         ordered = TRUE),
    Ferroptosis_score = as.numeric(Ferroptosis_score)
  ) %>%
  filter(!is.na(stage_event), !is.na(Ferroptosis_score))

#Box plot
pdf(file.path(output_dir, "FerroScores&TumorStage_box.pdf"),width = 9,height = 6)
grade_colors <- c("Stage I" = "#FFCCCC",  
                  "Stage II" = "#FFE6CC",   
                  "Stage III" = "#FFFFCC",    
                  "Stage IV" = "#CCFFCC",
                  "Stage X" = "#CCEEFF")  
ggplot(analysis_data_ss, aes(x = stage_event, y = Ferroptosis_score, fill = stage_event)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8 ) +
  geom_jitter(width = 0.15, alpha = 0.55, size = 3, shape = 21, color = "black",
              stroke = 0.55, aes(color = histologic_grade), show.legend = FALSE) +
  scale_fill_manual( values = grade_colors) +
  scale_color_manual(values = grade_colors)+
  labs(x = "Tumor Stage",
       y = "Ferroptosis Score") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        legend.position = "none"
  )
dev.off()

#Add median fitting curve
pdf(file.path(output_dir, "FerroScores&TumorStage_boxcurve.pdf"), width = 9, height = 6)
stage_levels <- levels(analysis_data_ss$stage_event)
median_data <- analysis_data_ss %>%
  group_by(stage_event) %>%
  summarise(median_score = median(Ferroptosis_score, na.rm = TRUE)) %>%
  mutate(stage_numeric = as.numeric(stage_event))
smooth_x <- seq(1, length(stage_levels), length.out = 100)
loess_fit <- loess(median_score ~ stage_numeric, data = median_data, span = 0.8)
smooth_y <- predict(loess_fit, newdata = data.frame(stage_numeric = smooth_x))
smooth_data <- data.frame(
  stage_numeric = smooth_x,
  median_score_smooth = smooth_y
)
x_breaks <- data.frame(
  stage_numeric = 1:length(stage_levels),
  stage_label = stage_levels
)
ggplot(analysis_data_ss, aes(x = stage_event, y = Ferroptosis_score, fill = stage_event)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.55, size = 3, shape = 21, color = "black",
              stroke = 0.55, show.legend = FALSE) +
  geom_point(data = median_data, 
             aes(x = stage_event, y = median_score), 
             shape = 18, size = 4, color = "#9370DB") +
  geom_line(data = smooth_data, 
            aes(x = stage_numeric, y = median_score_smooth), 
            color = "#9370DB", linewidth = 1.2, inherit.aes = FALSE) +
  scale_fill_manual(values = grade_colors) +
  scale_x_discrete(limits = stage_levels) +  
  labs(x = "Tumor Stage",
       y = "Ferroptosis Score") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        legend.position = "none")
dev.off()




#empirical curves for ferroptosis score & survival time across different stages
#Stage I Ferroptosis Score & Survival Time
FS_stageI <- read_csv("my_data/stage&Ferroscore&survivaltime_stageI.csv")

span_range <- seq(0.1, 1.0, by = 0.05)
gcv_values <- numeric(length(span_range))
trace_hat_values <- numeric(length(span_range))
for (i in seq_along(span_range)) {
  span_val <- span_range[i]
  loess_fit <- loess(survival_time ~ Ferroptosis_score, 
                     data = FS_stageI,
                     span = span_val,
                     degree = 2,
                     family = "gaussian")
  n <- nrow(FS_stageI)
  sse <- sum(loess_fit$residuals^2)
  trace_hat <- loess_fit$trace.hat
  gcv <- (sse / n) / ((1 - trace_hat / n)^2)
  gcv_values[i] <- gcv
  trace_hat_values[i] <- trace_hat
  cat(sprintf("%.2f\t%.4f\t\t%.2f\n", span_val, gcv, trace_hat))
}
optimal_idx <- which.min(gcv_values)
optimal_span <- span_range[optimal_idx]
min_gcv <- gcv_values[optimal_idx]
par(mar = c(5, 5, 4, 2) + 0.1)
plot(span_range, gcv_values, type = "b", pch = 16, col = "blue", lwd = 2,
     xlab = "Span Value", ylab = "GCV Value", cex.lab = 1.2, cex.axis = 1.1,
     main = "Optimal Span Selection Based on GCV Criterion", cex.main = 1.3)
abline(v = optimal_span, col = "red", lty = 2, lwd = 2)
points(optimal_span, min_gcv, col = "red", pch = 16, cex = 2)
text(optimal_span, min_gcv, 
     labels = sprintf("Optimal span: %.2f", optimal_span),
     pos = 3, col = "red", cex = 1.1)
grid()

pdf(file.path(output_dir, "StageI_FerroScores&SurvivalTime_empiricalcurve.pdf"), width = 9, height = 6)
FS_stageI$Ferroptosis_score <- as.numeric(FS_stageI$Ferroptosis_score)
FS_stageI$survival_time <- as.numeric(FS_stageI$survival_time)
p <- ggplot(FS_stageI, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18))
p + geom_smooth(method = "loess", span = 1, 
                se = FALSE, color = "#7291c0", fill = "#b8c8df", alpha = 0.2)
dev.off()

#Stage II Ferroptosis Score & Survival Time
FS_stageII <- read_csv("my_data/stage&Ferroscore&survivaltime_stageII.csv")

span_range <- seq(0.1, 1.0, by = 0.05)
gcv_values <- numeric(length(span_range))
trace_hat_values <- numeric(length(span_range))
for (i in seq_along(span_range)) {
  span_val <- span_range[i]
  loess_fit <- loess(survival_time ~ Ferroptosis_score, 
                     data = FS_stageII,
                     span = span_val,
                     degree = 2,
                     family = "gaussian")
  n <- nrow(FS_stageII)
  sse <- sum(loess_fit$residuals^2)
  trace_hat <- loess_fit$trace.hat
  gcv <- (sse / n) / ((1 - trace_hat / n)^2)
  gcv_values[i] <- gcv
  trace_hat_values[i] <- trace_hat
  cat(sprintf("%.2f\t%.4f\t\t%.2f\n", span_val, gcv, trace_hat))
}
optimal_idx <- which.min(gcv_values)
optimal_span <- span_range[optimal_idx]
min_gcv <- gcv_values[optimal_idx]
par(mar = c(5, 5, 4, 2) + 0.1)
plot(span_range, gcv_values, type = "b", pch = 16, col = "blue", lwd = 2,
     xlab = "Span Value", ylab = "GCV Value", cex.lab = 1.2, cex.axis = 1.1,
     main = "Optimal Span Selection Based on GCV Criterion", cex.main = 1.3)
abline(v = optimal_span, col = "red", lty = 2, lwd = 2)
points(optimal_span, min_gcv, col = "red", pch = 16, cex = 2)
text(optimal_span, min_gcv, 
     labels = sprintf("Optimal span: %.2f", optimal_span),
     pos = 3, col = "red", cex = 1.1)
grid()

pdf(file.path(output_dir, "StageII_FerroScores&SurvivalTime_empiricalcurve.pdf"), width = 9, height = 6)
FS_stageII$Ferroptosis_score <- as.numeric(FS_stageII$Ferroptosis_score)
FS_stageII$survival_time <- as.numeric(FS_stageII$survival_time)
p <- ggplot(FS_stageII, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18))
p + geom_smooth(method = "loess", span = 0.95, 
                se = FALSE, color = "#7291c0", fill = "#b8c8df", alpha = 0.2)
dev.off()

#Stage III&IV Ferroptosis Score & Survival Time
FS_stageIII <- read_csv("my_data/stage&Ferroscore&survivaltime_stageIII.csv")
FS_stageIV <- read_csv("my_data/stage&Ferroscore&survivaltime_stageIV.csv")
FS_stagelast <- rbind(FS_stageIII, FS_stageIV)

pdf(file.path(output_dir, "StageIII&IV_FerroScores&SurvivalTime_empiricalcurve.pdf"), width = 9, height = 6)
FS_stagelast$Ferroptosis_score <- as.numeric(FS_stagelast$Ferroptosis_score)
FS_stagelast$survival_time <- as.numeric(FS_stagelast$survival_time)
p <- ggplot(FS_stagelast, aes(x = Ferroptosis_score, y = survival_time)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  labs(x = "Ferroptosis Score", y = "Survival Time (day)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18))
p + geom_smooth(method = "loess",
                span = 1,
                se = FALSE,
                color = "#7291c0", 
                fill = "#b8c8df",
                alpha = 0.1)
dev.off()






