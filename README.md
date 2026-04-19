# NTHU AIS User Satisfaction Analysis

This repository contains the complete R analysis script for my IMBA thesis:
**"Evaluating User Satisfaction with NTHU's Academic Information System (AIS)"**

## Author
Dale John Baltazar  
IMBA in Technology Management  
National Tsing Hua University, Taiwan  

## Repository Contents

- `analysis.R` - Complete R script for all statistical analyses
- `LICENSE` - MIT License
- `README.md` - This file

## About This Research

This thesis examines user satisfaction with National Tsing Hua University's Academic Information System (AIS) using a comprehensive research model that includes:

- **Information Quality (IQ)**: Accuracy, completeness, relevance, and currency of information
- **System Quality (SQ)**: Reliability, accessibility, response time, and ease of use
- **Service Quality (SERVQ)**: Responsiveness, competence, and support
- **Mobile User Experience (MUX)**: Mobile interface design, navigation, and usability
- **Compatibility (CA)**: Integration with existing workflows and systems
- **User Satisfaction (US)**: Overall satisfaction and perceived value

## Reproducibility

All statistical analyses were performed in R and are **fully reproducible** using the provided `analysis.R` script.

### Requirements

- R (version 4.0 or higher recommended)
- R packages: `psych`, `dplyr`

### How to Run

1. Install required packages:
```r
install.packages(c("psych", "dplyr"))
```

2. Download your survey data CSV file (not included in this repository for privacy)

3. Update the file path in `analysis.R` (line 13):
```r
raw <- read.csv(
  "YOUR_PATH_HERE/Evaluating User Satisfaction with NTHU's Academic Information System.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
```

4. Run the entire script in R or RStudio

### What the Script Does

The `analysis.R` script performs:

1. **Data Loading & Cleaning**: Recodes Likert-scale responses to numeric values (1-5)
2. **Descriptive Statistics**: Mean, SD, skewness, kurtosis for all items
3. **Demographics**: Frequency tables for gender, age, education level, and usage frequency
4. **Reliability Analysis**: Cronbach's alpha for each construct
5. **Measurement Model Assessment**:
   - Composite Reliability (CR)
   - Average Variance Extracted (AVE)
   - Fornell-Larcker discriminant validity
   - HTMT discriminant validity
6. **Structural Model Assessment**:
   - Multiple regression with standardized coefficients
   - R², Adjusted R², F-statistic
   - Effect size (f²) for each predictor
   - Variance Inflation Factor (VIF) for multicollinearity
7. **Hypothesis Testing**: Support for H1-H5 (relationships between constructs and user satisfaction)

### Output

The script saves all results to `analysis_v2.RData`, which can be loaded in R for further analysis or visualization.

## Citation

If you use this code or methodology in your research, please cite:

```
Baltazar, D. J. (2026). Evaluating User Satisfaction with NTHU's Academic Information System (AIS).
IMBA Thesis, National Tsing Hua University, Taiwan.
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaboration, please open an issue in this repository.
