library(ggplot2)
library(dplyr)
library(sjPlot)
library(dplyr)

metabolomics_data = read.csv('/Users/angela/Downloads/sample_to_converter_hmp2.csv')

reads = read.csv("/Users/angela/Downloads/ismA_mapped_diagnosis.tsv", sep="\t")

diag <- unique(ibd_mgx[c(1,8)])

ibd_mgx <- unique(reads[,c(1,2,3,4)])

ibd_mgx_processed <- ibd_mgx %>%
  mutate(cpm = reads_mapped / total_reads * 1000000) 

total <- merge(metabolomics_data, ibd_mgx_processed, by="External.ID")

total <- total %>%
  mutate(Encoder = ifelse(cpm > 4, "Encoder", "Non-Encoder"))

df_summary <- total %>%
  group_by(Converter, Encoder) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Converter) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

# Create the stacked bar plot
ggplot(df_summary, aes(x = Converter, y = Count, fill = Encoder)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Encoder" = "purple", "Non-Encoder" = "yellow")) +  # Custom colors
  labs(title = "Percentage of Encoder vs. Non-Encoder by Converter Status",
       x = "Converter Status",
       y = "Percentage") +
  theme_minimal()

