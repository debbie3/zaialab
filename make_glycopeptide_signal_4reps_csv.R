# This script reads in replicate sets of glycresoft results and outputs
# a csv file with glycopeptide and total signal from each replicate.
# The output file is used in Jaccard and Tanimoto similarity scripts.

# Change input filepaths (10-13), total_signal column names (16-19), outfile names (63).

setwd("C:/Users/changd/Documents/Debbie/projects/similarity_R_script/20190418_make_glycan_signal_csv/5B8_SWZ_glycresoft_results/")

# read in glycresoft results
f1 <- read.csv("20190326_5B8eggL2_021119_tryp_HILIC_PRM_1ug_01_glycopeptides.csv", stringsAsFactors = FALSE)
f2 <- read.csv("20190326_5B8eggL2_021119_tryp_HILIC_PRM_1ug_02_glycopeptides.csv", stringsAsFactors = FALSE)
f3 <- read.csv("20190326_5B8eggL2_021119_tryp_HILIC_PRM_1ug_03_glycopeptides.csv", stringsAsFactors = FALSE)
f4 <- read.csv("20190326_5B8eggL2_021119_tryp_HILIC_PRM_1ug_04_glycopeptides.csv", stringsAsFactors = FALSE)

# rename total_signal column to reflect replicate name
names(f1)[names(f1) == "total_signal"] <- "5B8eggL2_tryp_PRM_01"
names(f2)[names(f2) == "total_signal"] <- "5B8eggL2_tryp_PRM_02"
names(f3)[names(f3) == "total_signal"] <- "5B8eggL2_tryp_PRM_03"
names(f4)[names(f4) == "total_signal"] <- "5B8eggL2_tryp_PRM_04"

# remove all rows except for those with protein_name == "P02763|A1AG1_HUMAN" for AGP
# or protein_name == "iav|ADJ41816|A/Phil/2-MA/1982_HA" for Phil82

# f1 <- f1[f1$protein_name == "P02763|A1AG1_HUMAN",]
# f2 <- f2[f2$protein_name == "P02763|A1AG1_HUMAN",]
# f3 <- f3[f3$protein_name == "P02763|A1AG1_HUMAN",]
# f4 <- f4[f4$protein_name == "P02763|A1AG1_HUMAN",]

f1 <- f1[f1$protein_name == "cus|SWZHA_5B8|A/Switzerland/9715293/2013",]
f2 <- f2[f2$protein_name == "cus|SWZHA_5B8|A/Switzerland/9715293/2013",]
f3 <- f3[f3$protein_name == "cus|SWZHA_5B8|A/Switzerland/9715293/2013",]
f4 <- f4[f4$protein_name == "cus|SWZHA_5B8|A/Switzerland/9715293/2013",]

# f1 <- f1[f1$protein_name == "cus|SWZHA|A/Switzerland/9715293/2013",]
# f2 <- f2[f2$protein_name == "cus|SWZHA|A/Switzerland/9715293/2013",]
# f3 <- f3[f3$protein_name == "cus|SWZHA|A/Switzerland/9715293/2013",]
# f4 <- f4[f4$protein_name == "cus|SWZHA|A/Switzerland/9715293/2013",]

# f1 <- f1[f1$protein_name == "iav|AFG99160|A/Phil/2/1982_HA",]
# f2 <- f2[f2$protein_name == "iav|AFG99160|A/Phil/2/1982_HA",]
# f3 <- f3[f3$protein_name == "iav|AFG99160|A/Phil/2/1982_HA",]
# f4 <- f4[f4$protein_name == "iav|AFG99160|A/Phil/2/1982_HA",]

# f1 <- f1[f1$protein_name == "iav|BRAZ78HA|A/Brazil/11/1978_HA",]
# f2 <- f2[f2$protein_name == "iav|BRAZ78HA|A/Brazil/11/1978_HA",]
# f3 <- f3[f3$protein_name == "iav|BRAZ78HA|A/Brazil/11/1978_HA",]
# f4 <- f4[f4$protein_name == "iav|BRAZ78HA|A/Brazil/11/1978_HA",]

# f1 <- f1[f1$protein_name == "iav|BRAZBS78|A/BS/Brazil/11/1978_HA",]
# f2 <- f2[f2$protein_name == "iav|BRAZBS78|A/BS/Brazil/11/1978_HA",]
# f3 <- f3[f3$protein_name == "iav|BRAZBS78|A/BS/Brazil/11/1978_HA",]
# f4 <- f4[f4$protein_name == "iav|BRAZBS78|A/BS/Brazil/11/1978_HA",]

# f1 <- f1[f1$protein_name == "iav|ACP41105|A/California/04/2009_HA",]
# f2 <- f2[f2$protein_name == "iav|ACP41105|A/California/04/2009_HA",]
# f3 <- f3[f3$protein_name == "iav|ACP41105|A/California/04/2009_HA",]
# f4 <- f4[f4$protein_name == "iav|ACP41105|A/California/04/2009_HA",]

# merge replicates
temp1 <- merge(f1, f2, by = "glycopeptide", all = TRUE)
temp2 <- merge(f3, f4, by = "glycopeptide", all = TRUE)

f <- merge(temp1, temp2, by = "glycopeptide", all = TRUE)
# f <- merge(temp1, f3, by = "glycopeptide", all = TRUE)



# sort all rows by peptide_start and neutral_mass
sort <- f[order(f$peptide_start.x.x, f$peptide_start.y.x, f$peptide_start.x.y,
                f$peptide_start.y.y, f$neutral_mass.x.x, f$neutral_mass.y.x, 
                f$neutral_mass.x.y, f$neutral_mass.y.y),]

# remove all columns except glycopeptide and total_signal
 vars <- names(f) %in% c("neutral_mass.x.x", "mass_accuracy.x.x", "ms1_score.x.x",
                        "ms2_score.x.x", "q_value.x.x", "start_time.x.x", "end_time.x.x",
                        "apex_time.x.x", "charge_states.x.x", "msms_count.x.x",
                        "peptide_start.x.x", "peptide_end.x.x", "protein_name.x.x",
                        "neutral_mass.y.x", "mass_accuracy.y.x", "ms1_score.y.x",
                        "ms2_score.y.x", "q_value.y.x", "start_time.y.x", "end_time.y.x",
                        "apex_time.y.x", "charge_states.y.x", "msms_count.y.x",
                        "peptide_start.y.x", "peptide_end.y.x", "protein_name.y.x",
                        "neutral_mass.x.y", "mass_accuracy.x.y", "ms1_score.x.y",
                        "ms2_score.x.y", "q_value.x.y", "start_time.x.y", "end_time.x.y",
                        "apex_time.x.y", "charge_states.x.y", "msms_count.x.y",
                        "peptide_start.x.y", "peptide_end.x.y", "protein_name.x.y",
                        "neutral_mass.y.y", "mass_accuracy.y.y", "ms1_score.y.y",
                        "ms2_score.y.y", "q_value.y.y", "start_time.y.y", "end_time.y.y",
                        "apex_time.y.y", "charge_states.y.y", "msms_count.y.y",
                        "peptide_start.y.y", "peptide_end.y.y", "protein_name.y.y")

# vars <- names(f) %in% c("neutral_mass.x", "mass_accuracy.x", "ms1_score.x", 
#                         "ms2_score.x", "q_value.x", "start_time.x", "end_time.x",
#                         "apex_time.x", "charge_states.x", "msms_count.x",
#                         "peptide_start.x", "peptide_end.x", "protein_name.x",
#                         "neutral_mass.y", "mass_accuracy.y", "ms1_score.y", 
#                         "ms2_score.y", "q_value.y", "start_time.y", "end_time.y",
#                         "apex_time.y", "charge_states.y", "msms_count.y",
#                         "peptide_start.y", "peptide_end.y", "protein_name.y",
#                         "neutral_mass", "mass_accuracy", "ms1_score", 
#                         "ms2_score", "q_value", "start_time", "end_time",
#                         "apex_time", "charge_states", "msms_count",
#                         "peptide_start", "peptide_end", "protein_name")


f <- f[!vars]

# replace NA with 0
f[is.na(f)] <- 0

# write output to a new file
dir.create("replicate_signal", showWarnings = FALSE)
write.csv(f, file = "replicate_signal/5B8eggL2_tryp_PRM_signal.csv", row.names = FALSE)

