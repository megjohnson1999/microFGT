#' Constants for microFGT example data generation
#'
#' This file contains constants used to generate realistic example data
#' representing female genital tract (FGT) microbiomes with different community 
#' state types (CSTs) and clinical parameters.

#' Taxonomic constants for community state types
#' Each CST has dominant, subdominant, and rare taxa
#' Format: list(taxon = c(King, Phylum, Class, Order, Family, Genus, Species), abundance = c(min, max))

#' @keywords internal
CST_I_TAXA <- list(
  # CST-I: Dominated by L. crispatus
  dominant = list(
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus crispatus"), 
         abundance = c(0.70, 0.95))
  ),
  subdominant = list(
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus jensenii"),
         abundance = c(0.05, 0.15)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus gasseri"),
         abundance = c(0.01, 0.10)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus iners"),
         abundance = c(0.01, 0.07))
  ),
  rare = list(
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Gardnerella", "Gardnerella vaginalis"),
         abundance = c(0.001, 0.03)),
    list(taxon = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Megasphaera", "Megasphaera sp."),
         abundance = c(0, 0.01)),
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Bifidobacterium", "Bifidobacterium sp."),
         abundance = c(0, 0.01)),
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella bivia"),
         abundance = c(0, 0.01)),
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella timonensis"),
         abundance = c(0, 0.01)),
    list(taxon = c("Bacteria", "Fusobacteria", "Fusobacteriia", "Fusobacteriales", "Leptotrichiaceae", "Sneathia", "Sneathia amnii"),
         abundance = c(0, 0.005))
  )
)

#' @keywords internal
CST_III_TAXA <- list(
  # CST-III: Dominated by L. iners
  dominant = list(
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus iners"), 
         abundance = c(0.60, 0.90))
  ),
  subdominant = list(
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Gardnerella", "Gardnerella vaginalis"),
         abundance = c(0.05, 0.20)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus crispatus"),
         abundance = c(0, 0.10)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus jensenii"),
         abundance = c(0, 0.10))
  ),
  rare = list(
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella bivia"),
         abundance = c(0, 0.05)),
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella timonensis"),
         abundance = c(0, 0.05)),
    list(taxon = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Megasphaera", "Megasphaera sp."),
         abundance = c(0, 0.03)),
    list(taxon = c("Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Shuttleworthia", "Shuttleworthia satelles"),
         abundance = c(0, 0.02)),
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Atopobium", "Atopobium vaginae"),
         abundance = c(0, 0.02)),
    list(taxon = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Dialister", "Dialister sp."),
         abundance = c(0, 0.02))
  )
)

#' @keywords internal
CST_IV_TAXA <- list(
  # CST-IV: Diverse community with anaerobes (BV-associated)
  dominant = list(
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Gardnerella", "Gardnerella vaginalis"), 
         abundance = c(0.15, 0.40))
  ),
  subdominant = list(
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella bivia"),
         abundance = c(0.05, 0.15)),
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella timonensis"),
         abundance = c(0.05, 0.15)),
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Atopobium", "Atopobium vaginae"),
         abundance = c(0.05, 0.15)),
    list(taxon = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Megasphaera", "Megasphaera sp."),
         abundance = c(0.05, 0.15)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus iners"),
         abundance = c(0, 0.20))
  ),
  rare = list(
    list(taxon = c("Bacteria", "Fusobacteria", "Fusobacteriia", "Fusobacteriales", "Leptotrichiaceae", "Sneathia", "Sneathia amnii"),
         abundance = c(0.01, 0.10)),
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Porphyromonadaceae", "Porphyromonas", "Porphyromonas asaccharolytica"),
         abundance = c(0.01, 0.10)),
    list(taxon = c("Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Shuttleworthia", "Shuttleworthia satelles"),
         abundance = c(0.01, 0.10)),
    list(taxon = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Dialister", "Dialister sp."),
         abundance = c(0.01, 0.10)),
    list(taxon = c("Bacteria", "Tenericutes", "Mollicutes", "Mycoplasmatales", "Mycoplasmataceae", "Mycoplasma", "Mycoplasma hominis"),
         abundance = c(0, 0.10)),
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Bifidobacterium", "Bifidobacterium sp."),
         abundance = c(0, 0.07)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Aerococcaceae", "Aerococcus", "Aerococcus christensenii"),
         abundance = c(0, 0.05)),
    list(taxon = c("Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Peptinophilaceae", "Peptinophilus", "Peptinophilus sp."),
         abundance = c(0, 0.05)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Streptococcaceae", "Streptococcus", "Streptococcus anginosus"),
         abundance = c(0, 0.05)),
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteriales", "Coriobacteriaceae", "Mobiluncus", "Mobiluncus mulieris"),
         abundance = c(0, 0.05))
  )
)

#' @keywords internal
CST_V_TAXA <- list(
  # CST-V: Dominated by L. jensenii
  dominant = list(
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus jensenii"), 
         abundance = c(0.70, 0.95))
  ),
  subdominant = list(
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus crispatus"),
         abundance = c(0.01, 0.15)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus gasseri"),
         abundance = c(0.01, 0.10)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus iners"),
         abundance = c(0.01, 0.07))
  ),
  rare = list(
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Gardnerella", "Gardnerella vaginalis"),
         abundance = c(0, 0.05)),
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella bivia"),
         abundance = c(0, 0.02)),
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Atopobium", "Atopobium vaginae"),
         abundance = c(0, 0.02)),
    list(taxon = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Megasphaera", "Megasphaera sp."),
         abundance = c(0, 0.01))
  )
)

#' @keywords internal
CST_DEFAULT <- list(
  # A mix of common taxa with no strong dominance
  dominant = list(
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus iners"), 
         abundance = c(0.20, 0.40)),
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Gardnerella", "Gardnerella vaginalis"), 
         abundance = c(0.10, 0.30)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus crispatus"), 
         abundance = c(0.05, 0.20))
  ),
  subdominant = list(
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella bivia"),
         abundance = c(0.05, 0.15)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus jensenii"),
         abundance = c(0.05, 0.15)),
    list(taxon = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus gasseri"),
         abundance = c(0.05, 0.15))
  ),
  rare = list(
    list(taxon = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Atopobium", "Atopobium vaginae"),
         abundance = c(0.01, 0.07)),
    list(taxon = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella timonensis"),
         abundance = c(0.01, 0.07)),
    list(taxon = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Megasphaera", "Megasphaera sp."),
         abundance = c(0.01, 0.07)),
    list(taxon = c("Bacteria", "Fusobacteria", "Fusobacteriia", "Fusobacteriales", "Leptotrichiaceae", "Sneathia", "Sneathia amnii"),
         abundance = c(0, 0.05)),
    list(taxon = c("Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Shuttleworthia", "Shuttleworthia satelles"),
         abundance = c(0, 0.05))
  )
)

#' Clinical parameter distributions by condition
#' Structure:
#'   parameter_name = list(
#'     condition = list(mean, sd, min, max)
#'   )
#' This allows generating realistic clinical parameters for different groups
#' @keywords internal
CLINICAL_PARAMS <- list(
  "pH" = list(
    "Healthy" = list(mean = 4.2, sd = 0.4, min = 3.8, max = 4.5),
    "BV" = list(mean = 5.2, sd = 0.5, min = 4.5, max = 6.0),
    "Intermediate" = list(mean = 4.7, sd = 0.3, min = 4.3, max = 5.0)
  ),
  "Nugent_Score" = list(
    "Healthy" = list(mean = 1.5, sd = 1.0, min = 0, max = 3),
    "BV" = list(mean = 7.5, sd = 1.5, min = 7, max = 10),
    "Intermediate" = list(mean = 5.0, sd = 1.0, min = 4, max = 6)
  ),
  "Amsel_Score" = list(
    "Healthy" = list(mean = 0.5, sd = 0.5, min = 0, max = 1),
    "BV" = list(mean = 3.0, sd = 0.8, min = 2, max = 4),
    "Intermediate" = list(mean = 1.5, sd = 0.5, min = 1, max = 2)
  ),
  "Shannon_Diversity" = list(
    "Healthy" = list(mean = 1.2, sd = 0.4, min = 0.6, max = 2.0),
    "BV" = list(mean = 2.4, sd = 0.3, min = 1.8, max = 3.0),
    "Intermediate" = list(mean = 1.8, sd = 0.3, min = 1.3, max = 2.3)
  ),
  "Simpson_Diversity" = list(
    "Healthy" = list(mean = 0.3, sd = 0.1, min = 0.1, max = 0.5),
    "BV" = list(mean = 0.7, sd = 0.1, min = 0.5, max = 0.9),
    "Intermediate" = list(mean = 0.5, sd = 0.1, min = 0.3, max = 0.7)
  )
)

#' CST to condition mapping
#' Maps community state types to likely clinical conditions
#' @keywords internal
CST_CONDITION_MAP <- list(
  "CST-I" = "Healthy",
  "CST-II" = "Healthy",
  "CST-III" = c("Healthy", "Intermediate"),
  "CST-IV" = c("BV", "Intermediate"),
  "CST-V" = "Healthy"
)

#' Taxonomic ranks
#' Standard taxonomic ranks used in microbiome data
#' @keywords internal
TAXONOMIC_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#' CST descriptions for documentation
#' @keywords internal
CST_DESCRIPTIONS <- list(
  "CST-I" = "Dominated by Lactobacillus crispatus",
  "CST-II" = "Dominated by Lactobacillus gasseri",
  "CST-III" = "Dominated by Lactobacillus iners",
  "CST-IV" = "Diverse anaerobic bacteria (associated with BV)",
  "CST-V" = "Dominated by Lactobacillus jensenii"
)