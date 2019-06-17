###################Get the raw data###########################
library(flowCore)
library(tidyverse)
library(stringr)

#May data
get_pro_data <- function() {
  
  #K1
  # pro1 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c01_sample-1_01_0_Day2_live.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c02_sample-1_01_0_Day4_live.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c03_sample-1_01_0_Day6_live.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c04_sample-1_01_0_Day8_live.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c05_sample-1_01_0_Day10_live.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c06_sample-1_01_0_Day11_live.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c16_sample-1_01_0_Day12_live.fcs", transformation=FALSE)
  # pro8 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c17_sample-1_01_0_Day14_live.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c18_sample-1_01_0_Day16_live.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c19_sample-1_01_0_Day18_live.fcs", transformation=FALSE)
  # pro11 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c20_sample-1_01_0_Day20_live.fcs", transformation=FALSE)
  
  #K2
  pro1 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c01_sample_01_0_Day0_live.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c02_sample_01_0_Day2_live.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c03_sample_01_0_Day4_live.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c04_sample_01_0_Day6_live.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c05_sample_01_0_Day8_live.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c07_sample_01_0_Day10_live.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c08_sample_01_0_Day11_live.fcs", transformation=FALSE)
  pro8 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c09_sample_01_0_Day12_live.fcs", transformation=FALSE)
  pro9 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c10_sample_01_0_Day14_live.fcs", transformation=FALSE)
  pro10 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c11_sample_01_0_Day16_live.fcs", transformation=FALSE)
  pro11 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c12_sample_01_0_Day18_live.fcs", transformation=FALSE)
  pro12 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c13_sample_01_0_Day20_live.fcs", transformation=FALSE)
  pro13 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c14_sample_01_0_Day22_live.fcs", transformation=FALSE)

  pro1_df = flowCore::exprs(pro1)
  pro2_df = flowCore::exprs(pro2)
  pro3_df = flowCore::exprs(pro3)
  pro4_df = flowCore::exprs(pro4)
  pro5_df = flowCore::exprs(pro5)
  pro6_df = flowCore::exprs(pro6)
  pro7_df = flowCore::exprs(pro7)
  pro8_df = flowCore::exprs(pro8)
  pro9_df = flowCore::exprs(pro9)
  pro10_df = flowCore::exprs(pro10)
  pro11_df = flowCore::exprs(pro11)
  pro12_df = flowCore::exprs(pro12)
  pro13_df = flowCore::exprs(pro13)
  
  g1 = data.frame(pro1_df) 
  g2 = data.frame(pro2_df) 
  g3 = data.frame(pro3_df) 
  g4 = data.frame(pro4_df) 
  g5 = data.frame(pro5_df) 
  g6 = data.frame(pro6_df) 
  g7 = data.frame(pro7_df) 
  g8 = data.frame(pro8_df) 
  g9 = data.frame(pro9_df) 
  g10 = data.frame(pro10_df) 
  g11 = data.frame(pro11_df)
  g12 = data.frame(pro12_df)
  g13 = data.frame(pro13_df)
  
  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f",
                 "CD90", "CD44", "CD41", "CD235ab")
  
  #other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B") #, "BCL11a")
  
  #cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  #DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194", "Event_Length")
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di",
              
              #"Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di") #, "Sm147Di")
  
  #"Nd148Di", "Nd150Di", "Er168Di", "I127Di",
  
  #"Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di", "Event_length")
  
  pro_marker = c(sur_marker, tran_factor)
  
  seq=c(1:26)
  
  day_levels = c(0,2,4,6,8,10,11,12,14,16,18,20,22) #c(2,4,6,8,10,11,12,14,16,18,20) # 
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    #K1
    # protein1 = bind_rows( g1[marker] %>% mutate(obstime = 2, col_id = 1:nrow(g1)),
    #                       g2[marker]%>% mutate(obstime = 4, col_id = 1:nrow(g2)),
    #                       g3[marker]%>% mutate(obstime = 6, col_id = 1:nrow(g3)),
    #                       g4[marker]%>% mutate(obstime = 8, col_id = 1:nrow(g4)),
    #                       g5[marker]%>% mutate(obstime = 10, col_id = 1:nrow(g5)),
    #                       g6[marker]%>% mutate(obstime = 11, col_id = 1:nrow(g6)),
    #                       g7[marker]%>% mutate(obstime = 12, col_id = 1:nrow(g7)),
    #                       g8[marker]%>% mutate(obstime = 14, col_id = 1:nrow(g8)),
    #                       g9[marker]%>% mutate(obstime = 16, col_id = 1:nrow(g9)),
    #                       g10[marker]%>% mutate(obstime = 18, col_id = 1:nrow(g10)),
    #                       g11[marker]%>% mutate(obstime = 20, col_id = 1:nrow(g11))) %>%
      
      #K2
    protein1 = bind_rows( g1[marker] %>% mutate(obstime = 0, col_id = 1:nrow(g1)),
                          g2[marker]%>% mutate(obstime = 2, col_id = 1:nrow(g2)),
                          g3[marker]%>% mutate(obstime = 4, col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(obstime = 6, col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(obstime = 8, col_id = 1:nrow(g5)),
                          g6[marker]%>% mutate(obstime = 10, col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(obstime = 11, col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(obstime = 12, col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(obstime = 14, col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(obstime = 16, col_id = 1:nrow(g10)),
                          g11[marker]%>% mutate(obstime = 18, col_id = 1:nrow(g11)),
                          g12[marker]%>% mutate(obstime = 20, col_id = 1:nrow(g12)),
                          g13[marker]%>% mutate(obstime = 22, col_id = 1:nrow(g13))) %>%

    mutate(obstime = factor(obstime, levels = day_levels))
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
  }   
  return(prot_list)
}

#Oct data
get_pro_data <- function() {
  
  #K1
  pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c01_K1_02_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c02_K1_02_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c03_K1_02_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c04_K1_02_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c05_K1_02_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c06_K1_02_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c07_K1_02_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c08_K1_02_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c09_K1_02_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c10_K1_02_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c15_K1_02_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c16_K1_02_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)
  # pro13 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c19_K1_02_0_0_MNCs_Ungated_Ungated.fcs", transformation=FALSE)
  # pro14 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c17_K1_02_0_0_Jurkat_Ungated_Ungated.fcs", transformation=FALSE)
  # 
  # # #K2
  # pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c01_K2_01_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c02_K2_01_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c03_K2_01_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c04_K2_01_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c05_K2_01_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c06_K2_01_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c07_K2_01_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  # pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c08_K2_01_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c09_K2_01_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c10_K2_01_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  # pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c15_K2_01_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  # pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c20_K2_01_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)
  # 
  #tibble(desc = as.character(pro13@parameters@data$desc), name = names(pro13@parameters@data$desc)) %>% separate(desc, c("Isotope", "Protein"), "_") %>% View()
  
  pro1_df = exprs(pro1)
  pro2_df = exprs(pro2)
  pro3_df = exprs(pro3)
  pro4_df = exprs(pro4)
  pro5_df = exprs(pro5)
  pro6_df = exprs(pro6)
  pro7_df = exprs(pro7)
  pro8_df = exprs(pro8)
  pro9_df = exprs(pro9)
  pro10_df = exprs(pro10)
  pro11_df = exprs(pro11)
  pro12_df = exprs(pro12)
  # pro13_df = exprs(pro13)
  # pro14_df = exprs(pro14)
  
  g1 = data.frame(pro1_df)
  g2 = data.frame(pro2_df)
  g3 = data.frame(pro3_df)
  g4 = data.frame(pro4_df)
  g5 = data.frame(pro5_df)
  g6 = data.frame(pro6_df)
  g7 = data.frame(pro7_df)
  g8 = data.frame(pro8_df)
  g9 = data.frame(pro9_df)
  g10 = data.frame(pro10_df)
  g11 = data.frame(pro11_df)
  g12 = data.frame(pro12_df)
  # g13 = data.frame(pro13_df)
  # g14 = data.frame(pro14_df)
  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", 
                 "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4")
  
  other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "MEF2C")
  
  cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")
  
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di", "Sm147Di", "Yb171Di",
              
              "Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di", "Nd146Di",
              
              "Nd148Di", "Nd150Di", "Er168Di", "I127Di",
              
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di") 
  
  pro_marker = c(sur_marker, other_marker, tran_factor, cell_cycle, DNA)
  
  seq=c(1:41)
  
  #day_levels = paste0("", c("MNCs", "Jurkats")) #paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20, "MNCs", "Jurkats"))
  day_levels = paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20))
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    protein1 = bind_rows( g1[marker] %>% mutate(Day = "0", col_id = 1:nrow(g1)), 
                          g2[marker]%>% mutate(Day = "2", col_id = 1:nrow(g2)),
                          g3[marker]%>% mutate(Day = "4", col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(Day = "6", col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(Day = "8", col_id = 1:nrow(g5)),
                          g6[marker]%>% mutate(Day = "10", col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(Day = "11", col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(Day = "12", col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(Day = "14", col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(Day = "16", col_id = 1:nrow(g10)), 
                          g11[marker]%>% mutate(Day = "18", col_id = 1:nrow(g11)),
                          g12[marker]%>% mutate(Day = "20", col_id = 1:nrow(g12))) %>% 
      # g13[marker]%>% mutate(Day = "MNCs", col_id = 1:nrow(g13)),
      # g14[marker]%>% mutate(Day = "Jurkats", col_id = 1:nrow(g14))) %>%
      mutate(Day = factor(Day, levels = day_levels))
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
    
  }  
  
  return(prot_list)
}

all_data = as_tibble(bind_rows(get_pro_data()))

all_data = all_data %>% mutate(cell = factor(paste(cell_id, obstime)))
all_data = all_data %>% mutate(capture = obstime)

gene_data0 = all_data %>% spread(cell.type, Intensity) 

#####gating: filter non-cells

mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  # reg = lm(data = gated_data, Ir191~Ir193)
  (predict(reg) + ir1)/2
}

gene_data = gene_data0 %>% filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
  filter(mean_Ir >= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5),
         mean_Ir <= sinh(8.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

gene_data %>% group_by(obstime) %>% summarise(count = n())

sur_marker = c("CD34", "CD71", "CD45RA", "CD123", "CD49f", 
               "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4",
               "CD36", "CD38")

other_marker = c("HBB", "HBA", "H3")

tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "TAL1", "GATA2",
                "RUNX1", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "Flag_Tag",
                "FLI1", "NFE2p45")

cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")

DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")

gated_data = gene_data0 %>%  #express_data0 %>%
  mutate(Day = obstime) %>% 
  select(-cell_id, -cell, -capture, -obstime) %>% #-cell_cycle, -HBB, -H3, -DNA, -mean_Ir, -mean_Pt
  mutate_if(is.numeric, ~ floor(.x)) %>% 
  filter((!Day %in% c("MNCs", "Jurkats")))

gated_data = gene_data %>%  #express_data0 %>%
  mutate(Day = obstime) %>% 
  select(-cell_id, -cell, -capture, -obstime, -cell_cycle, -HBB, -H3, -DNA, -mean_Ir, -mean_Pt) %>% #
  mutate_if(is.numeric, ~ floor(.x)) %>% 
  filter((!Day %in% c("MNCs", "Jurkats")))

prot_name = gated_data %>% select(-Day) %>% colnames() %>% as.character()

prot_list = data.frame(prot_name = prot_name, prot_index = 1:length(prot_name)) %>% as.tibble() %>% 
  mutate(prot_name = as.character(prot_name))

#summarise the gated_data
gated_data %>% group_by(Day) %>% summarise(count=n())

##################Load the cluster results#####################
#the data used for modelling
load(file = "Results_by_date/Dec2018/May/K3/sampled_cell_data.rda") #sub_data
load(file = "Results_by_date/Dec2018/May/K3/run_3.rda") #map_est4
#############extract the paprameters###############
  c = 25

  day_vals = sub_data$Day %>% unique()
  day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)
  
  prot_names = colnames(sub_data)
  prot_list = tibble(prot_name = prot_names, prot_index = 1:length(prot_names)) %>% filter(prot_name != "Day")
  
  run_fit = map_est4
  
  vals = run_fit$par 
  
  time = tibble(Day = day_index_0$day_index,  real_time= day_index_0$Day) 
  
  cluster = vals$theta %>% as.tibble() %>% pull(value) %>% as.tibble() %>% rename(vals = value) %>% mutate(cluster = 1:c)
  
  phi = vals$phi_s_cp %>% as.tibble() %>% mutate(cluster = 1:c) %>% gather(protein, std, -cluster) %>%
    mutate(protein = as.double(str_remove_all(protein, "V")))
  
  signal = vals$mu_s_cp %>% as.tibble() %>% mutate(cluster = 1:c) %>% gather(protein, vals, -cluster) %>%
    mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
    left_join(prot_list, by = c("protein" = "prot_index")) %>% 
    left_join(phi) 
  
  signal = signal %>% 
    group_by(prot_name) %>% 
    mutate(max_vals = max(vals), min_vals = min(vals), mid_vals = mean(vals), sd_vals = sd(vals)) %>% 
    ungroup()
  
  day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)
  
  mu_c_fit = vals$mu_s_cp %>% t()
  theta0_c_fit = vals$theta0_s_cp %>% t()
  phi_c_fit = vals$phi_s_cp %>% t()
  mu_bg_fit = vals$mu_bg %>% as.vector()

  st_fit = vals$scaling_t
  st_fit = c(1, st_fit)
  
  theta_fit = rep(cluster$vals, length(day_vals)) %>% 
    matrix(nrow = c)
  
  data_x = gated_data %>%
    as.tibble() %>% 
    mutate(Day = as.numeric(as.character(Day))) %>% 
    left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 
  
  data_x1 = data_x %>% select(-Day, -day_index) %>% as.matrix()
  
  day_indx = data_x %>% pull(day_index)
  
  all_cells = predict_cluster(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx)
  all_cells_prob = predict_cluster_prob(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx)

  ###########remove empty clusters#######
  print(table(all_cells)) %>% sum()
  
  #find the empty clusters
  low_cluster = table(all_cells) %>% as.tibble() %>% rename(cluster = all_cells) %>% 
    mutate(cluster = as.numeric(cluster)) %>% filter(n<=200) %>% pull(cluster)
  
  #find the missing clusters
  missing_cluster = base::setdiff(1:c, unique(all_cells) %>% as.tibble() %>% pull(V1))
  
  empty_cluster = c(low_cluster, missing_cluster)
  
  cluster_prob = all_cells_prob %>% as.matrix()
  
  cluster_prob = cluster_prob[, -(empty_cluster)]
  
  colnames(cluster_prob) <- paste0("p", 1:(c-length(empty_cluster)))
  
  old_cluster = 1:c
  
  old_cluster = old_cluster[-empty_cluster]
  
  new_cluster = data.frame(cluster = old_cluster, new_cluster = 1:length(old_cluster)) %>% as.tibble()
  
  order_cluster = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>%
    group_by(Day) %>% 
    mutate(prob = num/sum(num)) %>% ungroup() %>% 
    select(-num) %>% 
    group_by(cluster) %>%
    mutate(x = max(prob)) %>%
    filter(prob == x) %>%
    ungroup() %>% select(Day, cluster) %>% arrange(., Day) %>%
    select(cluster) %>% 
    filter(cluster %in% old_cluster) %>% 
    mutate(cluster = as.integer(cluster)) %>% 
    left_join(new_cluster) %>% 
    unique() %>% as.data.frame() %>% pull(cluster) %>% 
    as.integer()
  
  ordered_cluster = tibble(cluster = as.integer(order_cluster), state = 1:(c-length(empty_cluster)))
  
  library(ggrepel)
  ####################plot#########################
  p1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>% 
    filter(cluster %in% old_cluster) %>%
    left_join(new_cluster) %>%
    left_join(ordered_cluster) %>%
    #filter(state != 19, state != 20) %>% 
    ggplot(aes(x = prot_name, y= vals)) + 
    geom_point(aes(colour = prot_name)) +
    scale_y_log10()+
    #geom_text(aes(label = prot_name), size=3.0)+
    geom_text_repel(aes(label = prot_name), size=3.0) +
    facet_wrap(~state, nrow=2) +
    theme(legend.position = "none")
  
  p2 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
    mutate(prob = num/sum(num)) %>% ungroup() %>% 
    select(-num) %>% 
    spread(cluster, prob, fill = 0) %>% 
    gather(cluster, prob, -Day) %>% 
    mutate(cluster = as.numeric(cluster)) %>% 
    filter(cluster %in% old_cluster) %>% 
    left_join(new_cluster) %>% 
    left_join(ordered_cluster) %>% 
    ggplot(aes(x = Day, y = prob))+
    geom_point()+
    facet_wrap(~state, nrow=2)
  
  gridExtra::grid.arrange(p1, p2, ncol=1)
  
  gridExtra::grid.arrange(p1, p2, ncol=1)
  
  ###########plot the noise#############
  vals = map_est4$par 
  
  time = tibble(Day = day_index_0$day_index,  real_time= day_index_0$Day) 
  
  mu_bg_fit = vals$mu_bg %>% as.matrix() %>% t()
  
  st_fit = vals$scaling_t
  st_fit = c(1, st_fit) 
  
  scaler = st_fit %>% as.matrix()
  
  mu_n = scaler %*% mu_bg_fit  
  
  mu_n = mu_n %>% as.tibble() %>% mutate(Day = day_index_0$day_index) %>% left_join(time)
  colnames(mu_n) <- c(prot_list$prot_name, "Day", "obstime")
  
  mu_n %>% select(-obstime) %>% gather(protein, bg_noise, ATRX:TAL1) %>% 
    ggplot(aes(x=Day, y=bg_noise))+
    geom_point(aes(colour=protein))+
    facet_wrap(~protein, scales = "free_y")
  
  signal_mean = vals$mu_s_cp %>% as.matrix()
  signal_pi = vals$theta
  
  day_levels = time %>% pull(real_time)
  noise_predict = mu_n %>% select(-Day) %>% gather(prot_name, bg_noise, bCatenin:Thy1) %>% 
    mutate(bg_noise = asinh(bg_noise)) %>% 
    mutate(obstime = factor(obstime, levels = day_levels))
  ###################compare between runs#################
  prot_name = load(file = "Results_by_date/Dec2018/May/K3/sampled_cell_data.rda") %>% get() %>% select(-Day) %>% colnames()
  
  x = load(file = "Results_by_date/Dec2018/May/K3/run_1.rda") #map_est4
  data_1 = get(x)
  rm(x)
  
  x = load(file = "Results_by_date/Dec2018/May/K3/run_2.rda") #map_est4
  data_2 = get(x)
  rm(x)

  mu_1 = data_1$par$mu_s_cp
  colnames(mu_1) <- prot_name  
  
  mu_2 = data_2$par$mu_s_cp
  colnames(mu_2) <- prot_name  
  
  mu_2 %>% as.tibble() %>% mutate(run = 2, cluster = 1:n())
  mu_1 %>% as.tibble() %>% mutate(run = 1, cluster = 1:n()) %>% 
    gather(protein, Intensity1, ATRX:TAL1) %>% 
    left_join(mu_2 %>% as.tibble() %>% mutate(run = 2, cluster = 1:n())%>% 
                gather(protein, Intensity2, ATRX:TAL1), by = c("protein")) %>%
    pairs()
    
    ggplot(aes(x = Intensity1, y=Intensity2))+
    geom_point()+
    scale_y_log10()+
    scale_x_log10()+
    #geom_text(aes(label = protein), size=3.0)+
    facet_grid(cluster.x~cluster.y)
   ?pairs
    
  #########################find the connections################
  library(rstan)
  ######################find the nearest neighbours#############################
  corre_data = matrix(0, ncol = 5, nrow =ncol(cluster_prob)*ncol(cluster_prob))
    
  for(k in 1:ncol(cluster_prob)){ #!k %in% empty_cluster
      
    for(j in 1:ncol(cluster_prob)) {
        
      corr =0
      count_p1 = 0
      count_p2 = 0
        
      for(i in 1:nrow(cluster_prob)){
          
        p1 = cluster_prob[i, k]
        p2 = cluster_prob[i, j]
        pp = p1 + p2
          
        if(pp>=0.95) {
            
          if(p1 > 0.95) count_p1 = count_p1 + 1
          else if(p2 > 0.95) count_p2 = count_p2 + 1
          else corr = corr + 1
            
        }
      }
        
      index = ncol(cluster_prob) * (k-1) + j
        
      print(index)
        
      sum = count_p1 + count_p2 + corr
        
      corre_data[index, 1] =  k
      corre_data[index, 2] =  j
      corre_data[index, 3] =  count_p1/sum
      corre_data[index, 4] =  count_p2/sum
      corre_data[index, 5] =  corr/sum
    }
  }
    
knn_list = corre_data %>% as.tibble() %>% rename(cluster = V1, cluster_n = V2, p1 = V3, p2 = V4, p = V5) %>%
      select(-p1, -p2) %>% group_by(cluster) %>% arrange(cluster, desc(p)) %>% ungroup() %>% filter(p>=1e-3) %>% 
      select(cluster, cluster_n) %>%
      filter(cluster != cluster_n) %>% 
      left_join(new_cluster, by = c("cluster" = "new_cluster")) %>% 
      left_join(new_cluster, by = c("cluster_n" = "new_cluster")) %>% 
      select(-cluster.x, -cluster_n) %>% 
      rename(cluster = cluster.y, cluster_n = cluster.y.y) %>% 
      left_join(ordered_cluster, by = c("cluster" = "cluster")) %>% 
      left_join(ordered_cluster, by = c("cluster_n" = "cluster")) %>% 
      select(state.x, state.y)
  
    #bar plot for the density in each population  
  # corre_data %>% as.tibble() %>% rename(p1 = V1, p2 = V2, p = V3) %>%
  #   mutate(cluster = 2:(n()+1)) %>% gather(prob, count, p1:p) %>% 
  #   mutate(prob = factor(prob, levels = c("p1", "p2", "p"))) %>% 
  #   ggplot(aes(x = prob, y = count, fill = prob)) +
  #   geom_bar(stat="identity") +
  #   facet_wrap(~cluster, ncol = 6)
    
    # 
    # old_cluster = 1:c
    # 
    # old_cluster = old_cluster[-empty_cluster]
    # 
    # new_cluster = data.frame(cluster = old_cluster, new_cluster = 1:length(old_cluster)) %>% as.tibble()
    
    #plot the prob over time in each cluster, empty clusters are removed
    data_x %>% mutate(cluster = as.numeric(all_cells)) %>%
      left_join(new_cluster) %>% 
      group_by(Day, new_cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
      mutate(prob = num/sum(num)) %>% ungroup() %>% 
      select(-num) %>% 
      spread(new_cluster, prob, fill = 0) %>% 
      gather(new_cluster, prob, -Day) %>% 
      mutate(cluster = as.numeric(new_cluster)) %>% 
      left_join(ordered_cluster) %>% 
      ggplot(aes(x = Day, y = prob))+
      geom_point()+
      facet_wrap(~state, ncol=7)
    
    data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
      mutate(prob = num/sum(num)) %>% ungroup() %>% 
      select(-num) %>% 
      spread(cluster, prob, fill = 0) %>% 
      gather(cluster, prob, -Day) %>% 
      mutate(cluster = as.numeric(cluster)) %>% 
      filter(cluster %in% old_cluster) %>% 
      left_join(new_cluster) %>% 
      left_join(ordered_cluster) %>% 
      ggplot(aes(x = Day, y = prob))+
      geom_point()+
      facet_wrap(~state, nrow=2)
    
  ############dirichlet fit############  
    cmp_model = stan_model(file = "Software/cluster_model/dirich_stateV6.stan")
    
    data_t = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
      mutate(prob = num/sum(num)) %>% ungroup() %>% 
      select(-num) %>% 
      spread(Day, prob, fill = 0) %>%
      filter(cluster %in% old_cluster) %>% 
      left_join(new_cluster) %>% 
      left_join(ordered_cluster) %>% 
      arrange(state) %>% 
      #filter(!cluster %in% empty_cluster) %>% 
      select(`2`:`20`) %>% 
      as.matrix() %>% 
      t()
    
    n_state = ncol(data_t)
    n_time = nrow(data_t)
    n_rate = nrow(knn_list)
    
    time_points = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
      mutate(prob = num/sum(num)) %>% ungroup() %>% 
      select(-num) %>% 
      spread(Day, prob, fill = 0) %>%
      select(`2`:`20`) %>% colnames() %>% as.integer()
    
    fit_Yule_func = function(cmp_model, data_t, n_state, n_time, time_points, n_rate, knn_list){
      
      num_chains = 4
      
      data_for_stan = list(n_state = n_state, n_time = n_time, x_t_in = data_t, time_points = time_points, n_rate_in=n_rate, knn_weight = knn_list) 
      
      
      run_fit = sampling(cmp_model, iter = 600, warmup = 300, data = data_for_stan, chains = num_chains, thin = 1, cores = 14)
      
      run_fit
    }
    
    run_fit = fit_Yule_func(cmp_model, data_t, n_state, n_time, time_points, n_rate, knn_list)
    
    save(run_fit, file = "Results_by_date/Dec2018/May/K1/cell_death_model_2.rda")
    load(file = "Results_by_date/Dec2018/May/K1/cell_death_model_1_30.rda")
    traceplot(run_fit)
    
    plot(run_fit, pars = "lambda1")
    # 
    # png("Desktop/test.png", width = 1200, height = 1200)
    # test = rstan::extract(run_fit, pars = "alpha")
    # dim(test$alpha)
    # length(test)
    # dev.off()
    
    tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()

    alpha1 = rstan::summary(run_fit, par = "alpha1")
    alpha1$summary %>% View()

    alpha2 = rstan::summary(run_fit, par = "alpha2")
    # 
    # alpha$summary %>% as.tibble(rownames = "prob") %>%  select(prob, `50%`) %>% 
    #   mutate(prob = str_replace_all(prob, "alpha\\[|\\]", "")) %>% 
    #   separate(prob, c("To", "From"), ",") %>% 
    #   mutate(`50%` = if_else(To==From, 0, `50%`)) %>% #View()
    #   ggplot(aes(x = as.numeric(From), y = -as.numeric(To), fill= `50%`)) +
    #   geom_tile()
    #save(run_fit, file = "Desktop/Mar2018/fit_data_3.RData")
    #load(file = "Desktop/Mar2018/fit_data_1.RData")
    
    x_t_1 = rstan::summary(run_fit, par = "x_t_1") 
    
    x_t_1_new = x_t_1$summary
    
    predict = x_t_1_new %>% as.tibble(rownames = "prob") %>%  select(prob, `2.5%`, `50%`, `97.5%`)
    
    predict_data = predict %>% mutate(prob = str_replace_all(prob, "x_t_1\\[|\\]", "")) %>% 
      separate(prob, c("tp", "cluster"), ",") %>% mutate(tp = as.numeric(tp))
    
    time_index = data.frame(tp = c(1:length(time_points)), time = time_points) %>% 
      mutate(tp = as.numeric(tp)) %>% as.tibble()
    
    state_data = data_t %>% t() %>% as.tibble %>% mutate(cluster = 1:n()) %>% 
      gather(time, prob_val, -cluster) %>% mutate(time = as.numeric(time), cluster = as.character(cluster)) %>% 
      left_join(time_index) %>%  mutate(cluster = as.integer(cluster))
    
    ######order the cluster############
    ordered_cluster = tibble(cluster = order_cluster, state = 1:(c-length(empty_cluster)))
    
    order_cluster = state_data %>% 
      group_by(cluster) %>%
      mutate(x = max(prob_val)) %>%
      filter(prob_val == x) %>%
      ungroup() %>% select(time, cluster) %>% arrange(., time) %>%
      select(cluster) %>% unique() %>% as.data.frame() %>% pull(cluster) %>% 
      as.integer()
    
    #####################print the simulation data########################
    predict_data %>% 
      mutate(cluster = as.integer(cluster)) %>% 
      #left_join(state_data %>% mutate(new_clust = cluster, cluster = as.character(1:n()))) %>%
      #left_join(state_data) %>% 
      left_join(time_index) %>% 
      #left_join(ordered_cluster) %>% 
      #mutate(cluster = new_clust, tp = time) %>% 
      ggplot(aes(x = time, y = `50%`)) +
      geom_point(data = state_data %>% left_join(ordered_cluster), aes(x = time, y = prob_val, colour = "Experiment"), colour = "red") + 
      geom_line(aes(linetype = "Modelling"), size = 0.65) + 
      geom_ribbon(aes(ymin= `2.5%`, ymax =`97.5%`, size = "credible interval"), alpha = 0.2, fill ="blue") +
      #geom_point(data = state_data %>% mutate(cluster = as.numeric(as.character(cluster))) %>% gather(tp, prob_val, -cluster), 
      #            aes(x = as.numeric(tp), y = prob_val), colour = "red") + 
      #geom_point()+
      geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
      facet_wrap(~cluster, nrow = 2) + #, scales = "free_y"
      theme_classic() +
      theme(text = element_text(size=15)) +
      ylab("Proportion of Cells") +
      xlab("time") #+
     #theme(legend.title=element_blank())
    
    # all_data_1 %>%
    #   select(-Cell_id) %>% 
    #   select(Day, Ir191, Ir193, IdU, p_Rb) %>% 
    #   filter((!Day %in% c("Jurkats", "MNCs"))) %>% 
    #   as.tibble() %>% 
    #   mutate(Day = as.numeric(as.character(Day))) %>% 
    #   left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) %>% 
    #   mutate(cluster = all_cells) %>%
    #   filter(!cluster %in% empty_cluster) %>% 
    #   mutate(cluster = as.integer(cluster), cell_id = 1:n()) %>% 
    #   left_join(new_cluster) %>% select(-cluster) %>% 
    #   rename(cluster = new_cluster) %>% 
    #   left_join(ordered_cluster) %>% 
    #   select(-Day, -cluster, -day_index) %>%
    #   gather(Protein, Intensity,Ir191:p_Rb) %>% 
    #   mutate(protein = factor(Protein)) %>% 
    #   ggplot(aes(x=protein, y=asinh(Intensity), fill = protein)) +
    #   geom_violin(size=0.1, scale = "width", bw = 0.25) +
    #   geom_hline(yintercept = c(3.5), col = "red", lty = 2, alpha= 0.7) +
    #   geom_hline(yintercept = c(5.5), col = "purple", lty = 2, alpha= 0.7) +
    #   geom_hline(yintercept = c(6.15), col = "green4", lty = 2, alpha= 0.7) +
    #   facet_wrap(~state, nrow = 2) 
    
    
    signal %>% mutate(stds = vals*vals/(std*std), valss = (vals-min_vals)/(max_vals-min_vals)) %>% 
      mutate(cluster = as.integer(cluster)) %>% filter(!cluster %in% empty_cluster) %>% 
      left_join(new_cluster) %>% 
      select(-cluster) %>% rename(cluster = new_cluster) %>%
      left_join(ordered_cluster) %>% 
      ggplot(aes(x = prot_name, y= vals)) + 
      geom_point(aes(colour = prot_name)) +
      scale_y_log10()+
      #geom_errorbar(aes(ymin = asinh((vals-stds)/min), ymax = asinh((vals+stds)/min), colour = prot_name)) +
      geom_text(aes(label = prot_name), size=3.0)+
      facet_wrap(~state, ncol = 7) +
      theme(legend.position = "none")
    
    signal %>% mutate(stds = vals*vals/(std*std), valss = (vals-min_vals)/(max_vals-min_vals)) %>% 
      mutate(cluster = as.integer(cluster)) %>% filter(!cluster %in% empty_cluster) %>% 
      left_join(new_cluster) %>% 
      select(-cluster) %>% rename(cluster = new_cluster) %>% 
      #left_join(ordered_cluster) %>% 
      ggplot(aes(x = state, y= vals)) + 
      geom_point(aes(colour = state)) +
      scale_y_log10()+
      #geom_errorbar(aes(ymin = asinh((vals-stds)/min), ymax = asinh((vals+stds)/min), colour = prot_name)) +
      #geom_text(aes(label = prot_name), size=3.0)+
      facet_wrap(~prot_name, nrow = 2) +
      theme_classic() +
      theme(legend.position = "none")+
      xlab("cluster")
    
    prot_list
    
    xxxx = signal %>% mutate(stds = vals*vals/(std*std), valss = (vals-min_vals)/(max_vals-min_vals)) %>% 
      mutate(cluster = as.integer(cluster)) %>% filter(!cluster %in% empty_cluster) %>% 
      left_join(new_cluster) %>% 
      select(-cluster) %>% rename(cluster = new_cluster) %>% 
      left_join(ordered_cluster) %>% 
      select(state, protein, vals) %>% 
      left_join(prot_list, by = c("protein" = "prot_index")) %>% 
      select(-protein) %>% 
      spread(prot_name, vals) %>%
      select(-state) %>% 
      as.matrix() %>% asinh() 
    
    rownames(xxxx)<- c(1:20)
    xxxx%>% 
      pheatmap::pheatmap()
    
    ggplot(aes(x = state, y= vals)) + 
      geom_point(aes(colour = state)) +
      scale_y_log10()+
      #geom_errorbar(aes(ymin = asinh((vals-stds)/min), ymax = asinh((vals+stds)/min), colour = prot_name)) +
      #geom_text(aes(label = prot_name), size=3.0)+
      facet_wrap(~prot_name, nrow = 2) +
      theme_classic() +
      theme(legend.position = "none")+
      xlab("cluster")
    
    data_x %>% select(Day) %>% bind_cols(cluster_prob) %>% mutate(cell_id = 1:n()) %>% 
      gather(cluster, probability, p1:p14) %>% 
      select(cluster, probability) %>% 
      ggpairs(columns = 1, aes(colour = cluster, alpha = 0.4)) 
    
    ################knn network plot #######################
    library(RANN)
    library(tidygraph)
    library(ggraph)
    
    #######transform the lambda##########
    lambda_diff1 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
      mutate(prob = num/sum(num)) %>% ungroup() %>% 
      select(-num) %>% 
      spread(Day, prob, fill = 0) %>%
      filter(!cluster %in% empty_cluster) %>% 
      left_join(ordered_cluster) %>% 
      arrange(state) %>% 
      select(`2`:`11`) %>%
      mutate(cluster = 1:n()) %>% 
      gather(Day, lambda, -cluster) %>% 
      group_by(cluster) %>% 
      mutate(diff = max(lambda)-min(lambda)) %>% 
      ungroup() %>% select(cluster, diff) %>% 
      unique()
    
    lambda_diff2 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
      mutate(prob = num/sum(num)) %>% ungroup() %>% 
      select(-num) %>% 
      spread(Day, prob, fill = 0) %>%
      filter(!cluster %in% empty_cluster) %>% 
      left_join(ordered_cluster) %>% 
      arrange(state) %>% 
      select(`11`:`20`) %>%
      mutate(cluster = 1:n()) %>% 
      gather(Day, lambda, -cluster) %>% 
      group_by(cluster) %>% 
      mutate(diff = max(lambda)-min(lambda)) %>% 
      ungroup() %>% select(cluster, diff) %>% 
      unique()
    
    lambda_diff = lambda_diff1 %>% rename(diff1 = diff) %>% mutate(diff2 = lambda_diff2$diff)
    
    #find the transition rate matrix
    alpha_list1 = alpha1$summary %>% as.tibble(rownames = "prob") %>% select(prob, `50%`) %>% 
      mutate(prob = str_replace_all(prob, "alpha1\\[|\\]", "")) %>% 
      separate(prob, c("To", "From"), ",") %>% 
      mutate(`50%` = if_else(To==From, 0, `50%`)) %>% 
      #mutate(`2.5%` = if_else(To==From, 0, `2.5%`)) %>% 
      mutate(To = as.numeric(To), From = as.numeric(From)) %>% 
      rename(lambda1 = `50%`)
    
    alpha_list2 = alpha2$summary %>% as.tibble(rownames = "prob") %>% select(prob, `50%`) %>% 
      mutate(prob = str_replace_all(prob, "alpha2\\[|\\]", "")) %>% 
      separate(prob, c("To", "From"), ",") %>% 
      mutate(`50%` = if_else(To==From, 0, `50%`)) %>% 
      #mutate(`2.5%` = if_else(To==From, 0, `2.5%`)) %>% 
      mutate(To = as.numeric(To), From = as.numeric(From)) %>% 
      rename(lambda2 = `50%`)
    
    alpha_list = alpha_list1 %>% mutate(lambda2 = alpha_list2$lambda2) %>% 
      mutate(lambda = lambda1 + lambda2) %>% 
      select(To, From, lambda)
    
    alpha_list %>% View()
    
    connect_list = knn_list %>% rename(From = state.x, To = state.y)
    
    connect_list %>% 
      left_join(alpha_list) %>% 
      left_join(lambda_diff, by = c("From" = "cluster")) %>% 
      rename(diff_From = diff) %>% 
      left_join(lambda_diff, by = c("To" = "cluster")) %>% 
      rename(diff_To = diff) %>% 
      mutate(new_lambda = lambda * diff_From / (diff_To + 1e-10)) %>% 
      mutate(new_lambda = if_else(diff_To <= 1e-2 | diff_From <= 1e-2, 0, new_lambda)) %>% 
      select(-lambda, -diff_To, -diff_From) %>% 
      spread(To, new_lambda, fill = 0) %>% 
      View()
    
    tidy_knn = function(data_mat, k_val = 5){
      knn_network = nn2(data = data_mat, query = data_mat, k = k_val + 1)
      
      tbl_eucld   = knn_network$nn.dists %>%
        as_tibble() %>%
        magrittr::set_colnames(c("from", paste0("to", 1:k_val))) %>%
        mutate(from = 1:n()) %>%
        gather(to_num, eucl_dist, -from)
      
      tbl_edge    = knn_network$nn.idx %>%
        as_tibble() %>%
        magrittr::set_colnames(c("from", paste0("to", 1:k_val))) %>%
        gather(to_num, to, -from)
      
      tbl_edge %>%
        left_join(tbl_eucld, by = c("from", "to_num")) %>%
        select(-to_num) %>%
        as_tbl_graph() %>% left_join( data_mat %>% as_tibble() %>% mutate(id = as.character(1:n())), by=c('name'='id'))
    }
    
    # connect_list %>% 
    #   left_join(alpha_list) %>%
    #   rename(lambda = `50%`) %>% 
    #   filter(lambda >= 5e-2) %>% 
    #   as_tbl_graph() %>% 
    #   ggraph(layout = "fr") +
    #   geom_edge_link(aes(alpha = lambda), arrow = arrow(length = unit(6, 'mm'))) +
    #   geom_node_point(aes(col = name), size = 10) + geom_node_text(aes(label = name))
    # 
    ######### draw the full network ##############    
    network1 = connect_list %>% 
      left_join(alpha_list1) %>% 
      left_join(lambda_diff, by = c("From" = "cluster")) %>% 
      rename(diff_From1 = diff1) %>% 
      left_join(lambda_diff, by = c("To" = "cluster")) %>% 
      rename(diff_To1 = diff1) %>% 
      mutate(new_lambda1 = lambda1 * diff_From1 / (diff_To1+1e-10)) %>% 
      mutate(new_lambda1 = if_else(diff_To1 <= 1e-2 | diff_From1 <= 1e-2, 0, new_lambda1)) %>% #View()
      filter(new_lambda1 > 0.05) %>% select(From, To, new_lambda1) #%>% 
      # left_join(ordered_cluster, by = c("From" = "cluster")) %>% 
      # left_join(ordered_cluster, by = c("To" = "cluster")) %>% 
      # mutate(From = state.x, To = state.y) %>% 
      # select(-state.x, -state.y)
    
    network2 = connect_list %>%
      left_join(alpha_list2) %>%
      left_join(lambda_diff, by = c("From" = "cluster")) %>%
      rename(diff_From2 = diff2) %>%
      left_join(lambda_diff, by = c("To" = "cluster")) %>%
      rename(diff_To2 = diff2) %>%
      mutate(new_lambda2 = lambda2 * diff_From2 / (diff_To2+1e-10)) %>%
      mutate(new_lambda2 = if_else(diff_To2 <= 1e-2 | diff_From2 <= 1e-2, 0, new_lambda2)) %>%
      filter(new_lambda2 > 0.05) %>% select(From, To, new_lambda2) #%>% 
      #left_join(ordered_cluster, by = c("From" = "cluster")) %>% 
      #left_join(ordered_cluster, by = c("To" = "cluster")) %>% 
      #mutate(From = state.x, To = state.y) %>% 
      #select(-state.x, -state.y)
    
    connect_list %>% #left_join(ordered_cluster, by = c("From" = "cluster")) %>% 
      #left_join(ordered_cluster, by = c("To" = "cluster")) %>% 
      #mutate(From = state.x, To = state.y) %>% 
      #select(-state.x, -state.y) %>% 
      left_join(network1, by = c("From" = "From", "To" = "To")) %>% 
      left_join(network2, by = c("From" = "From", "To" = "To")) %>% 
      replace_na(replace = list(new_lambda1 = 0, new_lambda2 = 0)) %>% 
      mutate(new_lambda = (new_lambda1 + new_lambda2)/2) %>% 
      #mutate(new_lambda = if_else(new_lambda1 == 0, 0, new_lambda1 + new_lambda2)) %>% 
      #mutate(new_lambda = if_else(new_lambda2 == 0, 0, new_lambda1 + new_lambda2)) %>% 
      filter(new_lambda != 0) %>% 
      select(From, To, new_lambda, new_lambda1, new_lambda2) %>%
      # left_join(ordered_cluster, by = c("From" = "cluster")) %>% 
      # left_join(ordered_cluster, by = c("To" = "cluster")) %>% 
      # mutate(From = state.x, To = state.y) %>% 
      # select(-state.x, -state.y) %>% #pull(From) %>% unique()
      filter(!To %in% c(0), !From %in% c(0)) %>% 
      as_tbl_graph() %>% 
      activate(nodes) %>% 
      mutate(name = as.integer(name)) %>% 
      arrange(name) %>% 
      ggraph(layout = "auto") +#, node.positions = node_pos) +
      geom_edge_link(aes(alpha = new_lambda1), arrow = arrow(length = unit(5.5, 'mm')), colour = "red") +
      geom_edge_link(aes(alpha = new_lambda2), arrow = arrow(length = unit(5.5, 'mm')), colour = "blue") +
      geom_node_point(aes(col = factor(name, levels = 1:c)), size = 10) + geom_node_text(aes(label = name), size = 6.5) +
      theme_classic()+
      theme(legend.position = "none")+
      xlim(-5, 15)
    
    node_pos = matrix(c(1, 0, 10, 
                        2, 0, 5,
                        3, 0, 0,
                        4, 5, 0,
                        5, 5, 5,
                        6, 5, 10,
                        7, 10, -2.5,
                        8, 10, 2.5,
                        9, 10, 10,
                        10, 10, 5,
                        11, 15, 5,
                        12, 15, 10,
                        13, 17.5, 5,
                        14, 17.5, 10,
                        16, 20, 5,
                        17, 22.5, 5,
                        18, 25, 5), ncol = 3, byrow = T) %>% 
      as.tibble() %>% rename(name = V1, x = V2, y = V3) %>%
      #mutate(name = as.character(name)) %>% 
      arrange(name) %>% 
      select(x, y)
    ###############################################
    knn_list %>% as_tbl_graph() %>% 
      #activate(nodes) %>% 
      #left_join(node_pos) %>% 
      ggraph(layout = "auto") +
      geom_edge_link(arrow = arrow(length = unit(5, 'mm'))) +
      geom_node_point(aes(col = factor(name, levels = 1:c)), size = 8) + geom_node_text(aes(label = name))
    
    
    ## knn network all genes
    tdy_knn = tidy_knn(full_connect_list)
    
    #plot knn
    ggraph(tdy_knn, layout = "fr") +
      geom_edge_link(alpha = 0.1) +
      geom_node_point(aes(col = name, size = 0.1))
    
    connect_list %>% 
      left_join(alpha_list2) %>% 
      left_join(lambda_diff, by = c("From" = "cluster")) %>% 
      rename(diff_From = diff2) %>% 
      left_join(lambda_diff, by = c("To" = "cluster")) %>% 
      rename(diff_To = diff2) %>% 
      mutate(new_lambda = lambda2 * diff_From / (diff_To+1e-10)) %>% 
      mutate(new_lambda = if_else(diff_To <= 1e-2 | diff_From <= 1e-2, 0, new_lambda)) %>% #View()
      filter(new_lambda > 0.01) %>% 
      filter(!To %in% c(13, 20, 15), From != 13, From != 15) %>% 
      as_tbl_graph() %>% 
      #activate(nodes) %>% 
      #left_join(node_pos) %>% 
      ggraph(layout = "auto") +
      geom_edge_link(aes(alpha = new_lambda), arrow = arrow(length = unit(5, 'mm'))) +
      geom_node_point(aes(col = factor(name, levels = 1:c)), size = 8) + geom_node_text(aes(label = name))
    
  