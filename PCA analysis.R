# Install and load required packages
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggfortify)) install.packages("ggfortify")
if (!require(extrafont)) install.packages("extrafont")
library(ggplot2)
library(ggfortify)
library(extrafont)

# Register fonts for PDF output
font_import()  # Run this once to import system fonts (can take time)
loadfonts(device = "pdf")  # Register fonts for PDF device
# Verify available fonts
# print(fonts())  # Uncomment to check if "Times New Roman" is listed

# HC and PHC datasets
data_hc <- data.frame(
  Time = c(0.5, 1, 2, 4, 6),
  Reducing_Sugar = c(3085.1, 1351.11, 1822.19, 1180.72, 1220.81),
  Amino_Acid = c(1220.55, 931.8, 1285.1, 2209.09, 1149.21),
  X5.HMF = c(38.7, 22.74, 52.56, 25.64, 22.92),
  Furfuryl_Alcohol = c(51.29, 44.12, 58.79, 72.5, 11.8),
  Aldehyde = c(23.42, 31.26, 36.37, 28.79, 1.98),
  Ketone = c(2.81, 19.29, 4.08, 5.4, 8.12),
  Phenol = c(24.76, 9.25, 25.11, 31.52, 39.44),
  Amine = c(0, 0, 2.14, 3.12, 13.83),
  Thiophene = c(8.81, 5.34, 5.42, 6.89, 7.99),
  Furan = c(18.28, 2.89, 5.57, 2.18, 1.5),
  NCHC = c(18.37, 17.39, 16.63, 18.08, 24.81),
  pH = c(4.95, 4.99, 4.57, 4.64, 4.51),
  Humic_Acid = c(3.41, 5.66, 6.35, 6.76, 14.20)
)

data_phc <- data.frame(
  Time = c(0.5, 1, 2, 4, 6),
  Reducing_Sugar = c(408.94, 849.95, 829.91, 809.86, 609.4),
  Amino_Acid = c(3571.3, 3890.62, 4576.82, 4559.83, 2063.01),
  X5.HMF = c(20.28, 30.16, 32.81, 48.71, 12.42),
  Furfuryl_Alcohol = c(1.27, 24.76, 117.48, 226.15, 74.35),
  Aldehyde = c(21.23, 9.89, 3.42, 2.34, 3.05),
  Ketone = c(7.23, 6.17, 5.59, 7.61, 5.47),
  Phenol = c(18.59, 20.57, 27.35, 26.35, 23.76),
  Amine = c(24.5, 21.89, 17.59, 14.85, 28.64),
  Thiophene = c(6.07, 8.98, 17.79, 15.63, 11.7),
  Furan = c(0.32, 0.31, 0.57, 0.63, 0.48),
  NCHC = c(10.61, 14.61, 21.99, 27.29, 23.27),
  pH = c(8.87, 7.95, 8.13, 8.68, 7.61),
  Humic_Acid = c(7.97, 10.34, 11.90, 20.19, 14.67)
)

# Variables to analyze
vars_to_analyze <- c("Reducing_Sugar", "Amino_Acid", "X5.HMF", "Furfuryl_Alcohol", 
                     "Aldehyde", "Ketone", "Phenol", "Amine", "Thiophene", 
                     "Furan", "NCHC", "pH", "Humic_Acid")

# PCA for HC group
hc_data <- data_hc[, vars_to_analyze]
pca_hc <- prcomp(hc_data, scale. = TRUE)
summary(pca_hc)

# PCA for PHC group
phc_data <- data_phc[, vars_to_analyze]
pca_phc <- prcomp(phc_data, scale. = TRUE)
summary(pca_phc)

# Plot biplot for HC group
pdf("D:/Users/OneDrive/桌面/论文3/PCA/h1.pdf", width = 10, height = 8)
autoplot(pca_hc, data = hc_data, loadings = TRUE, loadings.label = TRUE, 
         loadings.label.size = 8, loadings.colour = '#bebada', 
         loadings.label.colour = "black", title = "PCA Biplot for HC Group (220°C, Hydrothermal)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, family = "Times"),
    axis.title = element_text(size = 24, family = "Times"),
    axis.text = element_text(size = 26, family = "Times")
  )
dev.off()

# Plot biplot for PHC group
pdf("D:/Users/OneDrive/桌面/论文3/PCA/ph4.pdf", width = 10, height = 8)
autoplot(pca_phc, data = phc_data, loadings = TRUE, loadings.label = TRUE, 
         loadings.label.size = 8, loadings.colour = '#fb8072', 
         loadings.label.colour = "black", title = "PCA Biplot for PHC Group (220°C, Hydrothermal)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, family = "Times"),
    axis.title = element_text(size = 24, family = "Times"),
    axis.text = element_text(size = 26, family = "Times")
  )
dev.off()




# Save PCA loadings and variance explained
write.csv(pca_hc$rotation, "D:/Users/OneDrive/桌面/论文3/PCA/pca_loadings_hc_220C.csv")
write.csv(pca_phc$rotation, "D:/Users/OneDrive/桌面/论文3/PCA/pca_loadings_phc_220C.csv")
write.csv(summary(pca_hc)$importance, "D:/Users/OneDrive/桌面/论文3/PCA/pca_variance_hc_220C.csv")
write.csv(summary(pca_phc)$importance, "D:/Users/OneDrive/桌面/论文3/PCA/pca_variance_phc_220C.csv")

