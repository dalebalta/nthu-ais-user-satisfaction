# ============================================================
# NTHU AIS User Satisfaction Analysis — IMBA Thesis
# Author: Dale John Baltazar
# Method: PLS-SEM via seminr (v2.4.2+)
# All analyses are fully reproducible using this script
# ============================================================

suppressPackageStartupMessages({
  library(psych)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(seminr)
})

# ─── SECTION 1: Load & Clean Data ────────────────────────────
raw <- read.csv(
  "Evaluating User Satisfaction with NTHU's Academic Information System.csv",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
cat("Total rows:", nrow(raw), "\n")
cat("Total cols:", ncol(raw), "\n")
colnames(raw)
head(raw, 3)

recode_likert <- function(x) {
  x2 <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  x2 <- trimws(tolower(x2))
  dplyr::case_when(
    x2 == "strongly agree"    ~ 5,
    x2 == "agree"             ~ 4,
    x2 == "neutral"           ~ 3,
    x2 == "disagree"          ~ 2,
    x2 == "strongly disagree" ~ 1,
    TRUE ~ NA_real_
  )
}

IQ    <- data.frame(lapply(raw[, 11:14], recode_likert))
SQ    <- data.frame(lapply(raw[, 15:18], recode_likert))
SERVQ <- data.frame(lapply(raw[, 19:21], recode_likert))
MUX   <- data.frame(lapply(raw[, 22:25], recode_likert))
CA    <- data.frame(lapply(raw[, 26:29], recode_likert))
US    <- data.frame(lapply(raw[, 30:33], recode_likert))

colnames(IQ)    <- paste0("IQ",    1:4)
colnames(SQ)    <- paste0("SQ",    1:4)
colnames(SERVQ) <- paste0("SERVQ", 1:3)
colnames(MUX)   <- paste0("MUX",   1:4)
colnames(CA)    <- paste0("CA",    1:4)
colnames(US)    <- paste0("US",    1:4)

cat("NA check — IQ:", sum(is.na(IQ)),
    " SQ:", sum(is.na(SQ)),
    " SERVQ:", sum(is.na(SERVQ)),
    " MUX:", sum(is.na(MUX)),
    " CA:", sum(is.na(CA)),
    " US:", sum(is.na(US)), "\n")

# ─── SECTION 2: Descriptive Statistics ───────────────────────
all_items <- cbind(IQ, SQ, SERVQ, MUX, CA, US)
desc <- psych::describe(all_items)
print(round(desc[, c("n","mean","sd","median","min","max","skew","kurtosis")], 3))

CS <- data.frame(
  IQ    = rowMeans(IQ,    na.rm = TRUE),
  SQ    = rowMeans(SQ,    na.rm = TRUE),
  SERVQ = rowMeans(SERVQ, na.rm = TRUE),
  MUX   = rowMeans(MUX,   na.rm = TRUE),
  CA    = rowMeans(CA,    na.rm = TRUE),
  US    = rowMeans(US,    na.rm = TRUE)
)
print(psych::describe(CS)[, c("n","mean","sd","min","max")])

# Demographics
gender_clean <- iconv(raw[, 8],  to = "ASCII//TRANSLIT", sub = "")
age_clean    <- iconv(raw[, 7],  to = "ASCII//TRANSLIT", sub = "")
lev_raw      <- iconv(raw[, 9],  to = "ASCII//TRANSLIT", sub = "")
freq_clean   <- iconv(raw[, 10], to = "ASCII//TRANSLIT", sub = "")

lev_simple <- ifelse(grepl("PhD",      lev_raw, ignore.case = TRUE), "PhD",
              ifelse(grepl("Master",   lev_raw, ignore.case = TRUE), "Master's",
              ifelse(grepl("ndergrad", lev_raw, ignore.case = TRUE), "Undergraduate",
              ifelse(grepl("Exchange", lev_raw, ignore.case = TRUE), "Exchange",
                     "Other"))))

cat("Gender:\n");         print(table(gender_clean))
cat("Age group:\n");      print(table(age_clean))
cat("Level of study:\n"); print(table(lev_simple))
cat("AIS frequency:\n");  print(table(freq_clean))

# ─── SECTION 3: Measurement Model — Alpha / CR / AVE ─────────
constructs_lst <- list(IQ=IQ, SQ=SQ, SERVQ=SERVQ, MUX=MUX, CA=CA, US=US)

alpha_vals <- sapply(constructs_lst, function(d)
  round(suppressWarnings(
    psych::alpha(d, check.keys=FALSE, warnings=FALSE)$total$raw_alpha), 3))
print(alpha_vals)

cr_ave_df <- data.frame(Construct=character(), Alpha=numeric(),
  CR=numeric(), AVE=numeric(), sqrt_AVE=numeric(), stringsAsFactors=FALSE)

for (nm in names(constructs_lst)) {
  fa_res <- suppressWarnings(tryCatch(
    psych::fa(constructs_lst[[nm]], nfactors=1, rotate="none", fm="minres", warnings=FALSE),
    error = function(e) psych::fa(constructs_lst[[nm]], nfactors=1, rotate="none", fm="pa")
  ))
  lam <- as.numeric(fa_res$loadings); lam[is.na(lam)] <- 0
  CR  <- round(sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2)), 3)
  AVE <- round(mean(lam^2), 3)
  cr_ave_df <- rbind(cr_ave_df, data.frame(
    Construct = nm, Alpha = alpha_vals[nm],
    CR = CR, AVE = AVE, sqrt_AVE = round(sqrt(AVE), 3),
    stringsAsFactors = FALSE))
}
print(cr_ave_df)

# Fornell-Larcker matrix
cor_mat <- cor(CS, use = "pairwise.complete.obs")
fl <- cor_mat
diag(fl) <- cr_ave_df$sqrt_AVE
cat("\nFornell-Larcker Matrix (diagonal = sqrt(AVE)):\n")
print(round(fl, 3))

# HTMT
get_htmt <- function(a_df, b_df) {
  a <- scale(a_df); b <- scale(b_df)
  cross_cors <- cor(a, b, use = "pairwise.complete.obs")
  within_a   <- cor(a, use = "pairwise.complete.obs")
  within_b   <- cor(b, use = "pairwise.complete.obs")
  mean_cross <- mean(abs(cross_cors), na.rm = TRUE)
  wa <- ifelse(ncol(a) > 1, mean(abs(within_a[lower.tri(within_a)])), 1)
  wb <- ifelse(ncol(b) > 1, mean(abs(within_b[lower.tri(within_b)])), 1)
  round(mean_cross / sqrt(wa * wb), 3)
}
clist  <- list(IQ=IQ, SQ=SQ, SERVQ=SERVQ, MUX=MUX, CA=CA, US=US)
cnames <- names(clist)
htmt_mat <- matrix(NA, 6, 6, dimnames = list(cnames, cnames))
for (i in 1:6) for (j in 1:6)
  if (i != j) htmt_mat[i, j] <- get_htmt(clist[[i]], clist[[j]])
cat("\nHTMT Matrix:\n"); print(round(htmt_mat, 3))

# ─── SECTION 4: PLS-SEM via seminr ───────────────────────────
data_items <- cbind(IQ, SQ, SERVQ, MUX, CA, US)

mm <- constructs(
  composite("IQ",    multi_items("IQ",    1:4)),
  composite("SQ",    multi_items("SQ",    1:4)),
  composite("SERVQ", multi_items("SERVQ", 1:3)),
  composite("MUX",   multi_items("MUX",   1:4)),
  composite("CA",    multi_items("CA",    1:4)),
  composite("US",    multi_items("US",    1:4))
)

sm_pls <- relationships(
  paths(from = c("IQ","SQ","SERVQ","MUX","CA"), to = "US")
)

fit <- estimate_pls(
  data              = data_items,
  measurement_model = mm,
  structural_model  = sm_pls
)

cat("\n=== PLS-SEM Summary ===\n")
summary_fit <- summary(fit)
print(summary_fit)

# Bootstrap inference (5000 resamples)
set.seed(2024)
boot <- bootstrap_model(fit, nboot = 5000, cores = 1, seed = 2024)
summary_boot <- summary(boot)
cat("\n=== Bootstrapped Path Coefficients ===\n")
print(summary_boot$bootstrapped_paths)

# R-squared (row = metric, col = construct)
r2_us  <- summary_fit$paths["R^2",    "US"]
r2_adj <- summary_fit$paths["AdjR^2", "US"]
cat("\nR-squared for US:",         round(r2_us,  3), "\n")
cat("Adjusted R-squared for US:", round(r2_adj, 3), "\n")

# ─── SECTION 5: Effect Sizes (f²) & VIF ──────────────────────
CS_std <- as.data.frame(scale(CS))
m_ols  <- lm(US ~ IQ + SQ + SERVQ + MUX + CA, data = CS_std)
sm_ols <- summary(m_ols)
r2full <- sm_ols$r.squared

preds <- c("IQ","SQ","SERVQ","MUX","CA")
cat(sprintf("\n%-6s %8s %-12s %8s\n", "Pred", "f2", "Magnitude", "VIF"))
for (p in preds) {
  others <- setdiff(preds, p)
  r2red  <- summary(lm(as.formula(paste("US~", paste(others, collapse="+"))),
                        data = CS_std))$r.squared
  f2     <- round((r2full - r2red) / (1 - r2full), 3)
  r2p    <- summary(lm(as.formula(paste(p, "~", paste(others, collapse="+"))),
                        data = CS_std))$r.squared
  vif    <- round(1 / (1 - r2p), 3)
  mag    <- ifelse(f2 >= 0.35, "Large",
            ifelse(f2 >= 0.15, "Medium",
            ifelse(f2 >= 0.02, "Small", "Negligible")))
  cat(sprintf("%-6s %8.3f %-12s %8.3f\n", p, f2, mag, vif))
}

# ─── SECTION 6: Hypothesis Summary Table (PLS-SEM) ───────────
bp <- as.data.frame(summary_boot$bootstrapped_paths)

hyp <- data.frame(
  Hypothesis = c("H1: IQ -> US","H2: SQ -> US","H3: SERVQ -> US",
                 "H4: MUX -> US","H5: CA -> US"),
  Beta       = round(bp[, "Original Est."],   3),
  Boot_Mean  = round(bp[, "Bootstrap Mean"],  3),
  SE         = round(bp[, "Bootstrap SD"],    3),
  T_Stat     = round(bp[, "T Stat."],         3),
  CI_Low     = round(bp[, "2.5% CI"],         3),
  CI_High    = round(bp[, "97.5% CI"],        3),
  p_value    = round(bp[, "Bootstrap P Val"], 4),
  Decision   = ifelse(bp[, "Bootstrap P Val"] < 0.05, "Supported", "Not Supported"),
  stringsAsFactors = FALSE
)
rownames(hyp) <- NULL
cat("\n=== Hypothesis Testing Results ===\n")
print(hyp)

# ─── SAVE ─────────────────────────────────────────────────────
save(raw, IQ, SQ, SERVQ, MUX, CA, US, all_items, CS,
     alpha_vals, cr_ave_df, cor_mat, fl, htmt_mat,
     m_ols, sm_ols, hyp, fit, boot, summary_fit, summary_boot,
     file = "NTHU_AIS_analysis_final.RData")
cat("\nAll results saved to NTHU_AIS_analysis_final.RData\n")
