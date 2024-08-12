library(tidyverse)
library(scales)
library(GenomicRanges)
library(VennDiagram)
library(geneRxCluster)

sites <- read_tsv('intSites.tsv')
sites$posid <- paste0(sites$chromosome, sites$strand, sites$position)

#---
s2m <- readr::read_tsv('sample2mass.tsv')
r <- bind_rows(lapply(split(sites, sites$internalSampleID), function(x){
  inputMass <- subset(s2m, GTSP == x$internalSampleID[1])$mass
  tibble(sample = x$internalSampleID[1],
         nSites = n_distinct(x$posid),
         Chao1 = unname(vegan::estimateR(x$estAbund)[2]),
         `inputMass (ng)` = inputMass,
         nSitesPerGenome = nSites / (inputMass/0.00646),
         Chao1PerGenome = Chao1 / (inputMass/0.00646))
}))

#---

sites.ig <- subset(sites, nearestFeatureDist == 0)
pairings <- read_tsv('pairings.tsv')

# Create Venn diagrams for myopaxon samples that come from the same iPSCs.

invisible(lapply(list(c('GTSP4888', 'GTSP4892', 'GTSP4893'), 
                      c('GTSP4891', 'GTSP4896', 'GTSP4897')), function(x){
   
  
   a <- subset(sites, internalSampleID == x[1])
   b <- subset(sites, internalSampleID == x[2])
   c <- subset(sites, internalSampleID == x[3])
                        
   w <- 5
   a$position1 <- a$position-w
   a$position2 <- a$position+w
   b$position1 <- b$position-w
   b$position2 <- b$position+w
   c$position1 <- c$position-w
   c$position2 <- c$position+w
                        
   a.g <- makeGRangesFromDataFrame(a, seqnames.field = 'chromosome', start.field = 'position1', end.field = 'position2', keep.extra.columns = TRUE)
   b.g <- makeGRangesFromDataFrame(b, seqnames.field = 'chromosome', start.field = 'position1', end.field = 'position2', keep.extra.columns = TRUE)
   c.g <- makeGRangesFromDataFrame(c, seqnames.field = 'chromosome', start.field = 'position1', end.field = 'position2', keep.extra.columns = TRUE)                 
       
   o1 <- findOverlaps(b.g, a.g)
   o1 <- o1[! duplicated(queryHits(o1))]
   o1 <- o1[! duplicated(subjectHits(o1))]
   o1 <- data.frame(o1)
   b.g[o1$queryHits]$posid <- a.g[o1$subjectHits]$posid
   
   o2 <- findOverlaps(c.g, a.g)
   o2 <- o2[! duplicated(queryHits(o2))]
   o2 <- o2[! duplicated(subjectHits(o2))]
   o2 <- data.frame(o2)
   c.g[o2$queryHits]$posid <- a.g[o2$subjectHits]$posid
       
  
   # Consider sites not found in iPCSs and see if they can be merged 
   # with each other to increase overlap between non-iPCS samples.
   b.g1 <- b.g[o1$queryHits]
   c.g1 <- c.g[o2$queryHits]        
   b.g2 <- b.g[-o1$queryHits]
   c.g2 <- c.g[-o2$queryHits]
   
   o3 <- findOverlaps(b.g2, c.g2)
   o3 <- o3[! duplicated(queryHits(o3))]
   o3 <- o3[! duplicated(subjectHits(o3))]
   o3 <- data.frame(o3)
   
   b.g2[o3$queryHits]$posid <- c.g2[o3$subjectHits]$posid
   
   s <- list()
   s[[a$subject[1]]] <- unique(a.g$posid)
   s[[b$subject[1]]] <- unique(c(b.g1$posid, b.g2$posid))
   s[[c$subject[1]]] <- unique(c(c.g1$posid, c.g2$posid))
   names(s) <- sub('^p', '', names(s))
             
   venn.diagram(s, fill = c("gold", "green", "red"), 
                   alpha = c(0.5, 0.5, 0.5), paste0(a$subject[1], '3ways_venn.png'), 
                   units = 'in', height = 5, width = 5, imagetype = 'png',
                   resolution = 300, cat.just=list(c(-2, 1.5), c(1.5,4), c(0.5,-3.5)))
                   file.remove(list.files(pattern = 'log'))
  }))



sharedAbunds <- bind_rows(lapply(split(pairings, 1:nrow(pairings)), function(x){
  
  IPSC <- sub('^p', '', unique(subset(sites, internalSampleID == x$IPSC)$subject))
  Myopaxon <- sub('^p', '', unique(subset(sites, internalSampleID == x$Myopaxon)$subject))
  
  a <- subset(sites, internalSampleID == x$IPSC)
  b <- subset(sites, internalSampleID == x$Myopaxon)
  
  w <- 5
  a$position1 <- a$position-w
  a$position2 <- a$position+w
  b$position1 <- b$position-w
  b$position2 <- b$position+w
  
  a.g <- makeGRangesFromDataFrame(a, seqnames.field = 'chromosome', start.field = 'position1', end.field = 'position2', keep.extra.columns = TRUE)
  b.g <- makeGRangesFromDataFrame(b, seqnames.field = 'chromosome', start.field = 'position1', end.field = 'position2', keep.extra.columns = TRUE)
  
  o <- findOverlaps(b.g, a.g)
  o <- o[! duplicated(queryHits(o))]
  o <- o[! duplicated(subjectHits(o))]
  o <- data.frame(o)
  
  # Determine shared clonal abundance.
  b.hits <- subset(b, posid %in% b.g[o$queryHits]$posid)
  a.hits <- subset(a, posid %in% a.g[o$subjectHits]$posid)
  
  sharedAbunds <- tibble(iPSC = IPSC, 
                         iPSC_percentCellsShared =  sprintf("%.2f%%", (sum(a.hits$estAbund) / sum(a$estAbund))*100),
                         Myopaxon = Myopaxon, 
                         Myopaxon_percentCellsShared =  sprintf("%.2f%%", (sum(b.hits$estAbund) / sum(b$estAbund))*100))
  
  # Update myopaxon posids with IPSC ids.
  b.g[o$queryHits]$posid <- a.g[o$subjectHits]$posid
  
  # Create a venn diagram.
  s <- list(IPSC = a.g$posid, Myopaxon = b.g$posid)
  names(s) <- c(IPSC, Myopaxon)
  venn.diagram(s, fill = c("lightblue", "green"), 
               alpha = c(0.5, 0.5), paste0(IPSC, '_', Myopaxon, '_venn.png'), 
               units = 'in', height = 5, width = 5, imagetype = 'png',
               resolution = 300, cat.just=list(c(-2,0), c(1,-5)))
  file.remove(list.files(pattern = 'log'))
  
  
  # Calculate nearest gene frequencies (10KB).
  a <- subset(a, abs(nearestFeatureDist) <= 10000)
  b <- subset(b, abs(nearestFeatureDist) <= 10000)
  
  p <- bind_rows(lapply(unique(c(a$nearestFeature, b$nearestFeature)), function(x){
         tibble(g = x, 
                a = n_distinct(subset(a, nearestFeature == x)$posid)/n_distinct(a$posid), 
                b = n_distinct(subset(b, nearestFeature == x)$posid)/n_distinct(b$posid))
       }))
  
  p$d <- p$a - p$b / sqrt(2)
  
  lab1 <- arrange(p, d)[1:5,]$g
  lab2 <- arrange(p, desc(d))[1:5,]$g
  
  p <- bind_rows(subset(p, ! g %in% c(lab1, lab2)), subset(p, g %in% lab1), subset(p, g %in% lab2))
  p$z <- ifelse(p$g %in% c(lab1, lab2), p$g, 'Other')
  p$z <- factor(p$z, levels = c(lab1, lab2, 'Other'))
  
  reds <- c('#770101', '#9e3e36', '#c26f68', '#e1a09d', '#fed3d3')
  blues <- c('#1a01b2', '#4c3cc8', '#7167db', '#9691ed', '#bcbbfc')
  
  set.seed(1)
  o <- ggplot(p, aes(a, b, fill = z)) + 
       theme_bw()+
       labs(x = 'iPSC gene frequency', y = 'Myopaxon gene frequency') +
       scale_x_continuous(limits = c(-0.001, max(c(p$a, p$b))+0.0005)) +
       scale_y_continuous(limits = c(-0.001, max(c(p$a, p$b))+0.0005))+
       scale_fill_manual(name = '', values = c(blues, reds, 'gray50')) +
       geom_jitter(width = 0.0001, height = 0.0001, size = 3.5, color = 'black', shape = 21, alpha = 0.9) +
       geom_abline(intercept = 0, slope = 1, color = 'red') +
       theme(panel.grid.minor = element_blank(),
             axis.text = element_text(size = 12),
             axis.title=element_text(size=14),
             legend.key.size = unit(2, 'line'),
             legend.text=element_text(size=14),
             axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
             axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
       guides(fill = guide_legend(override.aes = list(size = 7))) +
       ggtitle(paste(IPSC, ' / ', Myopaxon))
  
  ggsave(o, file = paste0(IPSC, '_', Myopaxon, '_geneFreq.png'), units = 'in', 
         width = 8, height = 5, dpi = 300)
  
  n <- 16
  colors <- c('gray90', grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n))
  
  b <- arrange(b, desc(relAbund))
  d <- tibble(g = b[1:n,]$nearestFeature,
              r = b[1:n,]$relAbund)
  d <- bind_rows(d, tibble(g = 'LowAbund', r = 100-sum(d$r)))
  d$g <- factor(d$g, levels = c('LowAbund', b[1:n,]$nearestFeature))
  d$s <- Myopaxon
  
  p <- ggplot(d, aes(s, r/100, fill = g)) +
       theme_bw() +
       scale_fill_manual(name = 'Gene', values = colors) +
       geom_col() +
       scale_y_continuous(expand=c(0,0), labels = scales::percent_format(accuracy = 1)) +
       scale_x_discrete(expand=c(0,0)) +
       labs(x = '', y = '') +
       theme(axis.text = element_text(size = 14))
  
  ggsave(p, file = paste0(IPSC, '_', Myopaxon, '_relAbubnd.png'), units = 'in', 
         height = 5, width = 4, dpi = 300)
  
  return(sharedAbunds)
}))



saveRDS(sharedAbunds, file = 'sharedAbunds.rds')

#----

sites.ig <- group_by(sites.ig, internalSampleID) %>% mutate(sampleSitesNum = n_distinct(posid)) %>% ungroup()

d <- bind_rows(lapply(split(sites.ig, paste(sites.ig$subject, sites.ig$nearestFeature)), function(x){
       tibble(subject = x$subject[1], gene = x$nearestFeature[1], geneIntFreq = n_distinct(x$posid) / x$sampleSitesNum[1])
     })) 

o <- data.frame(spread(d, gene, geneIntFreq))
rownames(o) <- o$subject
o <- o[,c(2:length(o))]
o[is.na(o)] <- 0
pca <- prcomp(o, center = TRUE, scale. = TRUE)
library(ggbiplot)
ggbiplot(pca)



o <- as.data.frame(t(o))
names(o) <- o[1,]
o <- o[2:nrow(o),]
o[is.na(o)] <- 0
pca <- prcomp(o, center = TRUE, scale. = TRUE)

#----

topGenes <- unique(unlist(lapply(split(pairings, 1:nrow(pairings)), function(x){
  IPSC <- sub('^p', '', unique(subset(sites, internalSampleID == x$IPSC)$subject))
  Myopaxon <- sub('^p', '', unique(subset(sites, internalSampleID == x$Myopaxon)$subject))
  
  a <- subset(sites, internalSampleID == x$IPSC)
  b <- subset(sites, internalSampleID == x$Myopaxon)
  b <- subset(b, abs(b$nearestFeatureDist) <= 10000)
  i <- n_distinct(b$posid)
  
  d <- group_by(b, nearestFeature) %>%
       summarise(f = n_distinct(posid) / i) %>%
       ungroup() %>% arrange(desc(f))
  
  # Includes ties.
  slice_max(d, f, n = 10)$nearestFeature
})))



d <- bind_rows(lapply(split(pairings, 1:nrow(pairings)), function(x){
  IPSC <- sub('^p', '', unique(subset(sites, internalSampleID == x$IPSC)$subject))
  Myopaxon <- sub('^p', '', unique(subset(sites, internalSampleID == x$Myopaxon)$subject))
  
  a <- subset(sites, internalSampleID == x$IPSC)
  b <- subset(sites, internalSampleID == x$Myopaxon)
  b <- subset(b, abs(b$nearestFeatureDist) <= 10000)
  i <- n_distinct(b$posid)
  
  bind_rows(lapply(topGenes, function(g){
    tibble(s = Myopaxon, g = g, f = n_distinct(subset(b, nearestFeature == g)$posid) / i)
  }))
}))

o <- data.frame(sort(table(subset(d, f > 0)$g), decreasing = TRUE))
d <- left_join(d, o, by = c('g'='Var1')) 
d <- arrange(d, desc(Freq), desc(f))
d$g <- factor(d$g, levels = rev(unique(d$g)))
d$s <- factor(d$s, levels = gtools::mixedsort(unique(d$s)))

d[d$f == 0,]$f <- NA

p <- ggplot(d, aes(s, g, fill = f)) + 
     geom_tile(color = 'black') +
     scale_fill_gradient2(name = 'Frequency',
                          low = 'green2', mid = 'yellow', high = 'red', 
                          na.value = 'white', midpoint = 0.0025) +
    labs(x = '', y = '') +
    scale_x_discrete(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0))

ggsave(p, file = paste0('top_genes.png'), units = 'in', 
       height = 15, width = 8, dpi = 300)





samples <- readr::read_tsv('published_gt_sample_data/published_gt_samples.tsv')
samples <- samples[grepl('lenti', samples$Vector, ignore.case = TRUE),]

z <- subset(samples, timePointDays <= 31)
z[z$CellType %in% c('BM CD34+', 'CD34', 'CD34_0214', 'CD34_1014', 'CD34+ Cells'),]$CellType <- 'CD34+'
z[z$CellType %in% c('T cells', 'Tcells:CAR+', 'Tcells:CAR+CD8-', 'Tcells:CAR+CD8+'),]$CellType <- 'Tcells'
z[z$CellType %in% c('Whole blood', 'Whole Blood'),]$CellType <- 'Whole blood'
z[z$CellType %in% c('BM', 'Bone Marrow'),]$CellType <- 'Bone marrow'

# Reduce the number of trials to make plotting reasonable.
z <- z[! z$Trial %in% c('PSMA', 'Rivella', 'UPENN_CART-UPCC32816', 'UPENN_CART19_CLL', 'Gill_CART19_18415'),]

z[grepl('beta', z$Trial, ignore.case = TRUE),]$Trial <- 'betaThalassemia'
z[grepl('WAS', z$Trial, ignore.case = TRUE),]$Trial <- 'WAS'
z[grepl('UPENN_CART19_ALL', z$Trial, ignore.case = TRUE),]$Trial <- 'CART19 ALL'

z <- z[! grepl('hit', z$Timepoint),]
z$Timepoint <- toupper(z$Timepoint)

pubSites <- read_tsv('published_gt_sample_data/published_gt_sites.tsv')
pubSites$posid <- paste0(pubSites$chromosome, pubSites$strand, pubSites$position)

pubSites <- group_by(pubSites, internalSampleID) %>%
         mutate(sampleSites = n_distinct(posid)) %>%
         ungroup() %>%
         filter(sampleSites >= 100 & internalSampleID %in% z$GTSP)

d <- bind_rows(lapply(split(pubSites, pubSites$internalSampleID), function(x){
       tibble(GTSP = x$internalSampleID[1], 
              percentOnco = n_distinct(subset(x, abs(x$nearestOncoFeatureDist) < 50000)$posid) / n_distinct(x$posid))

     }))

z <- left_join(z, d, by = 'GTSP')
z <- z[! is.na(z$percentOnco),]

# Rename patients.
n <- 1
z <- bind_rows(lapply(split(z, z$Trial), function(x){
       bind_rows(lapply(split(x, x$Patient), function(x2){
         x2$Patient <- paste0('p', n)
         n <<- n + 1
         x2
        }))
     }))


# Read in project sites and calculate % near oncogenes.
sites <- read_tsv('intSites.tsv')
sites$posid <- paste0(sites$chromosome, sites$strand, sites$position)

d2 <- bind_rows(lapply(split(sites, sites$internalSampleID), function(x){
  tibble(Trial = 'Myopaxon', Patient = x$subject[1], Timepoint = x$timePoint[1], 
         CellType = x$cellType[1],GTSP = x$internalSampleID[1], 
         percentOnco = n_distinct(subset(x, abs(x$nearestOncoFeatureDist) < 50000)$posid) / n_distinct(x$posid))
}))

z <- z[,names(d2)]


# Merge project sites and published sites.
pd <- bind_rows(z, d2)

# Remove samples which were sequenced more than once, use latest sequencing results.
pd <- pd[! duplicated(pd$GTSP),]
pd <- pd[! pd$GTSP %in% c('GTSP0644', 'GTSP0652', 'GTSP0647', 'GTSP0654', 'GTSP0630', 'GTSP0886'),]

pd <- group_by(pd, Trial) %>% 
      arrange(percentOnco) %>% 
      ungroup() %>% 
      arrange(factor(Trial, levels = c("Myopaxon", "CGD", "betaThalassemia", "WAS",  "CART19 ALL")))

pd$Patient <- sub('^pMyo', 'Myo', pd$Patient)
pd$Patient <- sub('^pIPSC', 'IPSC', pd$Patient)
pd$Patient <- sub('^Myopaxon', 'Myopax', pd$Patient)
pd$label <- paste0(pd$Patient, '-', pd$Timepoint)
pd$label <- factor(pd$label, levels = unique(pd$label))
pd$Trial <- factor(pd$Trial, levels = unique(pd$Trial))

p <- ggplot(pd, aes(label, percentOnco, fill = CellType)) +
     theme_bw() +
     scale_fill_manual(name = 'Cell type', 
                       values = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(7)) +
     scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
     geom_col(position = position_dodge2(width = 1, preserve = "single")) +
     labs(x = '', y = 'Percent integration near oncogenes') +
     theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
           axis.text.y = element_text(size = 10), 
           axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
           legend.position="bottom",
           panel.grid.major.x = element_blank(),
           panel.grid.major.y = element_line( size=.1, color="black" ),
           strip.background = element_rect(fill='white'),
           strip.text = element_text(color = 'black')) + 
     facet_wrap(Trial~., scales = 'free_x', ncol = 2)

ggsave(p, file = 'pNearOnco.png', dpi = 300,units = 'in', height = 8, width = 12)



m <- list()
m[['iPSC vs Myopaxon']] =  
  matrix(c(n_distinct(subset(sites[  grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites[  grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(sites[! grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites[! grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) <= 50000)$posid)),
          byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('iPSC', 'Myopaxon')))

m[['iPSC vs published trials']] = 
  matrix(c(n_distinct(subset(sites[grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites[grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(pubSites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(pubSites, abs(nearestOncoFeatureDist) <= 50000)$posid)),
         byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('myo', 'pub')))

m[['Myopaxon vs published trials']] = 
  matrix(c(n_distinct(subset(sites[! grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites[! grepl('IPSC', sites$subject),], abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(pubSites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(pubSites, abs(nearestOncoFeatureDist) <= 50000)$posid)),
         byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('myo', 'pub')))

m[['(iPSC + Myopaxon) vs published trials']] = 
  matrix(c(n_distinct(subset(sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites, abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(pubSites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(pubSites, abs(nearestOncoFeatureDist) <= 50000)$posid)),
         byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('myo', 'pub')))

WAS_sites <- subset(pubSites, internalSampleID %in% subset(z, Trial == 'WAS')$GTSP)
m[['(iPSC + Myopaxon) vs WAS']] = 
  matrix(c(n_distinct(subset(sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites, abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(WAS_sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(WAS_sites, abs(nearestOncoFeatureDist) <= 50000)$posid)),
         byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('myo', 'pub')))

CGD_sites <- subset(pubSites, internalSampleID %in% subset(z, Trial == 'CGD')$GTSP)
m[['(iPSC + Myopaxon) vs CGD']] = 
  matrix(c(n_distinct(subset(sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites, abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(CGD_sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(CGD_sites, abs(nearestOncoFeatureDist) <= 50000)$posid)),
         byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('myo', 'pub')))

betaThalassemia_sites <- subset(pubSites, internalSampleID %in% subset(z, Trial == 'betaThalassemia')$GTSP)
m[['(iPSC + Myopaxon) vs betaThalassemia']] = 
  matrix(c(n_distinct(subset(sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites, abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(betaThalassemia_sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(betaThalassemia_sites, abs(nearestOncoFeatureDist) <= 50000)$posid)),
         byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('myo', 'pub')))

CARTsites <- subset(pubSites, internalSampleID %in% subset(z, Trial == 'CART19 ALL')$GTSP)
m[['(iPSC + Myopaxon) vs CART19 ALL']] = 
  matrix(c(n_distinct(subset(sites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(sites, abs(nearestOncoFeatureDist) <= 50000)$posid),
           n_distinct(subset(CARTsites, abs(nearestOncoFeatureDist) >  50000)$posid),
           n_distinct(subset(CARTsites, abs(nearestOncoFeatureDist) <= 50000)$posid)),
         byrow = FALSE, ncol = 2, dimnames = list(c('Not near', 'Near'), c('myo', 'pub')))

md <- bind_rows(mapply(function(x, n){
  tibble(test = n,  
         A_nearOnco = sprintf("%.2f%%", (x[2,1]/sum(x[,1]))*100), 
         B_nearOnco = sprintf("%.2f%%", (x[2,2]/sum(x[,2]))*100),
         pVal = fisher.test(x)$p.value)
}, m, names(m), SIMPLIFY = FALSE))

saveRDS(md, file = 'FisherTests.rds')


calc_uc50 <- function(x){
  stopifnot(is.vector(x) & is.numeric(x))
  x <- x[order(x)]
  accum <- sapply(1:length(x), function(i){sum(x[1:i])})
  length(accum[accum >= sum(x)/2])
}



sampleTable <- bind_rows(lapply(split(sites, sites$internalSampleID), function(x){
                 tibble('Sample ID' = x$internalSampleID[1], 'Sample Name' = sub('^p', '', x$subject[1]), 
                        'Cell Type' = x$cellType[1], Reads = sum(x$reads), 'Inferred Cells' = sum(x$estAbund),
                        'Unique Sites' = n_distinct(x$posid), 'S. Chao1' = unname(vegan::estimateR(x$estAbund)[2]),
                       'Shannon' = vegan::diversity(x$estAbund, index = 'shannon'), 'Gini' = ineq::ineq(x$estAbund, type="Gini"),
                       'Simpson' = vegan::diversity(x$estAbund, index = 'simpson'), 'UC50' = calc_uc50(x$estAbund))
                }))

saveRDS(sampleTable, file = 'sampleTable.rds')


sampleTable2 <- readr::read_tsv('sampleTable.tsv')




