# NTHU AIS User Satisfaction Analysis - IMBA Thesis
# Author: Dale John Baltazar
# All analyses are fully reproducible using this script

suppressPackageStartupMessages({
  library(psych)
  library(dplyr)
})

# Read CSV dataset (update path to your local file)
raw <- read.csv(
  "Evaluating User Satisfaction with NTHU's Academic Information System.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Total rows:", nrow(raw), "\n")
cat("Total cols:", ncol(raw), "\n")
colnames(raw)
head(raw, 3)

# Data Cleaning and Recoding
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

colnames(IQ)    <- paste0("IQ", 1:4)
colnames(SQ)    <- paste0("SQ", 1:4)
colnames(SERVQ) <- paste0("SERVQ", 1:3)
colnames(MUX)   <- paste0("MUX", 1:4)
colnames(CA)    <- paste0("CA", 1:4)
colnames(US)    <- paste0("US", 1:4)

cat("NA counts — IQ:",  sum(is.na(IQ)),
    " SQ:",    sum(is.na(SQ)),
    " SERVQ:", sum(is.na(SERVQ)),
    " MUX:",   sum(is.na(MUX)),
    " CA:",    sum(is.na(CA)),
    " US:",    sum(is.na(US)), "\n")

CS <- data.frame(
  IQ    = rowMeans(IQ,    na.rm = TRUE),
  SQ    = rowMeans(SQ,    na.rm = TRUE),
  SERVQ = rowMeans(SERVQ, na.rm = TRUE),
  MUX   = rowMeans(MUX,   na.rm = TRUE),
  CA    = rowMeans(CA,    na.rm = TRUE),
  US    = rowMeans(US,    na.rm = TRUE)
)

describe(CS)[, c("n", "mean", "sd", "min", "max")]

# Descriptive Statistics
all_items <- cbind(IQ, SQ, SERVQ, MUX, CA, US)
desc <- describe(all_items)
print(round(desc[, c("n", "mean", "sd", "median", "min", "max", "skew", "kurtosis")], 3))

gender_clean <- iconv(raw[, 8],  to = "ASCII//TRANSLIT", sub = "")
age_clean    <- iconv(raw[, 7],  to = "ASCII//TRANSLIT", sub = "")
lev_raw      <- iconv(raw[, 9],  to = "ASCII//TRANSLIT", sub = "")
freq_clean   <- iconv(raw[, 10], to = "ASCII//TRANSLIT", sub = "")

lev_simple <- ifelse(grepl("PhD", lev_raw, ignore.case = TRUE), "PhD",
  ifelse(grepl("Master", lev_raw, ignore.case = TRUE), "Master's",
    ifelse(grepl("ndergrad", lev_raw, ignore.case = TRUE), "Undergraduate",
      ifelse(grepl("Exchange", lev_raw, ignore.case = TRUE), "Exchange", "Other"))))

cat("Gender:\n"); print(table(gender_clean))
cat("Age group:\n"); print(table(age_clean))
cat("Level of study:\n"); print(table(lev_simple))
cat("AIS frequency:\n"); print(table(freq_clean))

# Measurement Model
constructs <- list(IQ = IQ, SQ = SQ, SERVQ = SERVQ, MUX = MUX, CA = CA, US = US)
alpha_vals <- sapply(constructs, function(d)
  round(suppressWarnings(psych::alpha(d, check.keys = FALSE, warnings = FALSE)$total$raw_alpha), 3))
print(alpha_vals)

cr_ave_df <- data.frame(Construct = character(), Alpha = numeric(), CR = numeric(),
  AVE = numeric(), sqrt_AVE = numeric(), stringsAsFactors = FALSE)
fa_loadings <- list()

for (nm in names(constructs)) {
  fa_res <- suppressWarnings(tryCatch(
    psych::fa(constructs[[nm]], nfactors = 1, rotate = "none", fm = "minres", warnings = FALSE),
    error = function(e) psych::fa(constructs[[nm]], nfactors = 1, rotate = "none", fm = "pa")))
  lam <- as.numeric(fa_res$loadings)
  lam[is.na(lam)] <- 0
  fa_loadings[[nm]] <- round(lam, 3)
  CR  <- round(sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2)), 3)
  AVE <- round(mean(lam^2), 3)
  cr_ave_df <- rbind(cr_ave_df, data.frame(Construct = nm, Alpha = alpha_vals[nm],
    CR = CR, AVE = AVE, sqrt_AVE = round(sqrt(AVE), 3), stringsAsFactors = FALSE))
}
print(cr_ave_df)

cor_mat <- cor(CS, use = "pairwise.complete.obs")
fl <- cor_mat
diag(fl) <- cr_ave_df$sqrt_AVE
cat("Fornell-Larcker Matrix:\n")
print(round(fl, 3))

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

clist <- list(IQ = IQ, SQ = SQ, SERVQ = SERVQ, MUX = MUX, CA = CA, US = US)
htmt_mat <- matrix(NA, 6, 6, dimnames = list(names(clist), names(clist)))
for (i in 1:6) for (j in 1:6) {
  if (i != j) htmt_mat[i, j] <- get_htmt(clist[[i]], clist[[j]])
}
print(round(htmt_mat, 3))

# Structural Model
m_std <- lm(scale(US) ~ scale(IQ) + scale(SQ) + scale(SERVQ) + scale(MUX) + scale(CA), data = CS)
sm <- summary(m_std)

cat("R-squared:", round(sm$r.squared, 3), "\n")
cat("Adjusted R-squared:", round(sm$adj.r.squared, 3), "\n")
fstat <- sm$fstatistic
cat("F(", fstat[2], ",", fstat[3], ") =", round(fstat[1], 3),
  " p =", round(pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE), 4), "\n\n")
print(round(sm$coefficients, 4))

preds <- c("IQ", "SQ", "SERVQ", "MUX", "CA")
r2full <- sm$r.squared
for (p in preds) {
  others <- setdiff(preds, p)
  r2red <- summary(lm(as.formula(paste("scale(US) ~",
    paste(paste0("scale(", others, ")"), collapse = "+"))), data = CS))$r.squared
  f2 <- round((r2full - r2red) / (1 - r2full), 3)
  mag <- ifelse(f2 >= 0.35, "Large", ifelse(f2 >= 0.15, "Medium",
    ifelse(f2 >= 0.02, "Small", "Negligible")))
  r2p <- summary(lm(as.formula(paste(p, "~", paste(others, collapse = "+"))), data = CS))$r.squared
  vif <- round(1 / (1 - r2p), 3)
  cat(p, "| f2:", f2, mag, "| VIF:", vif, "\n")
}

keys <- c("scale(IQ)", "scale(SQ)", "scale(SERVQ)", "scale(MUX)", "scale(CA)")
coefs <- sm$coefficients[keys, ]
hyp <- data.frame(
  Hypothesis = c("H1: IQ->US", "H2: SQ->US", "H3: SERVQ->US", "H4: MUX->US", "H5: CA->US"),
  Beta = round(coefs[, 1], 3), SE = round(coefs[, 2], 3),
  t = round(coefs[, 3], 3), p = round(coefs[, 4], 4),
  Decision = ifelse(coefs[, 4] < 0.05, "Supported", "Not Supported"),
  stringsAsFactors = FALSE
)
rownames(hyp) <- NULL
print(hyp)

save(raw, IQ, SQ, SERVQ, MUX, CA, US, all_items, CS, alpha_vals, fa_loadings,
  cr_ave_df, cor_mat, fl, htmt_mat, m_std, sm, hyp, file = "analysis_v2.RData")
