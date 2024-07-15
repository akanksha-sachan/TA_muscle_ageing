# FastIIB
FastIIB <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/FastIIB.csv", row.names = 1, header = TRUE)
FastIIB[1:4]
dim(FastIIB)
subset_FastIIB <- subset(FastIIB, p_val < 0.01 & abs(avg_log2FC) > 1)
dim(subset_FastIIB) # 322  5

cat(rownames(subset_FastIIB), sep = "\n")
writeLines(rownames(subset_FastIIB), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/FastIIB.txt")

# FastIIX
FastIIX <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/FastIIX.csv", row.names = 1, header = TRUE)
FastIIX[1:4]
dim(FastIIX)
subset_FastIIX <- subset(FastIIX, p_val < 0.01 & abs(avg_log2FC) > 1)
dim(subset_FastIIX) # 154  5

cat(rownames(subset_FastIIX), sep = "\n")
writeLines(rownames(subset_FastIIX), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/FastIIX.txt")

# FAPs
FAPs <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/FAPs.csv", row.names = 1, header = TRUE)
FAPs[1:4]
dim(FAPs)
subset_FAPs <- subset(FAPs, p_val < 0.01 & abs(avg_log2FC) > 1)
dim(subset_FAPs) # 383   5

cat(rownames(subset_FAPs), sep = "\n")
writeLines(rownames(subset_FAPs), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/FAPs.txt")

# Skeleton MuSc
Skeleton_MuSc <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Skeleton_MuSc.csv", row.names = 1, header = TRUE)
Skeleton_MuSc[1:4]
dim(Skeleton_MuSc)
subset_Skeleton_MuSc <- subset(Skeleton_MuSc, p_val < 0.01 & abs(avg_log2FC) > 1)
dim(subset_Skeleton_MuSc) # 200   5

cat(rownames(subset_Skeleton_MuSc), sep = "\n")
writeLines(rownames(subset_Skeleton_MuSc), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Skeleton_MuSc.txt")

# MTJ
MTJ <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/MTJ.csv", row.names = 1, header = TRUE)
MTJ[1:4]
dim(MTJ)
subset_MTJ <- subset(MTJ, p_val < 0.01 & avg_log2FC > 1)
dim(subset_MTJ) # 87   5

cat(rownames(subset_MTJ), sep = "\n")
writeLines(rownames(subset_MTJ), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/MTJ.txt")

# Pericyte
Pericyte <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Pericyte.csv", row.names = 1, header = TRUE)
Pericyte[1:4]
dim(Pericyte)
subset_Pericyte <- subset(Pericyte, p_val < 0.01 & avg_log2FC > 1)
dim(subset_Pericyte) # 181   5

cat(rownames(subset_Pericyte), sep = "\n")
writeLines(rownames(subset_Pericyte), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Pericyte.txt")

# EC
EC <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/EC.csv", row.names = 1, header = TRUE)
EC[1:4]
dim(EC)
subset_EC <- subset(EC, p_val < 0.01 & avg_log2FC > 1)
dim(subset_EC) # 104   5

cat(rownames(subset_EC), sep = "\n")
writeLines(rownames(subset_EC), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/EC.txt")

# Macrophages
Macrophages <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Macrophages.csv", row.names = 1, header = TRUE)
Macrophages[1:4]
dim(Macrophages)
subset_Macrophages <- subset(Macrophages, p_val < 0.01 & avg_log2FC > 1)
dim(subset_Macrophages) # 53   5

cat(rownames(subset_Macrophages), sep = "\n")
writeLines(rownames(subset_Macrophages), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Macrophages.txt")

# NMJ
NMJ <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/NMJ.csv", row.names = 1, header = TRUE)
NMJ[1:4]
dim(NMJ)
subset_NMJ <- subset(NMJ, p_val < 0.01 & avg_log2FC > 1)
dim(subset_NMJ) # 65   5

cat(rownames(subset_NMJ), sep = "\n")
writeLines(rownames(subset_NMJ), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/NMJ.txt")

# Tendon
Tendon <- read.csv("/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Tendon.csv", row.names = 1, header = TRUE)
Tendon[1:4]
dim(Tendon)
subset_Tendon <- subset(Tendon, p_val < 0.01 & avg_log2FC > 1)
dim(subset_Tendon) # 108   5

cat(rownames(subset_Tendon), sep = "\n")
writeLines(rownames(subset_Tendon), "/ix/djishnu/Zarifeh/ML_MM/Aditi/result/ClusterSpecific_Cells/DEGs_ALL/Tendon.txt")
