Pathway ZIKV
================
Joan M. Valls Cuevas
5/23/2019

Read in the data and relabel the columns

``` r
my_dat <- data.frame(read.delim("hDC_AR_rlog.txt"))
my_dat2 <- data.frame(my_dat[,9:17], row.names = my_dat[,1])
colnames(my_dat2) <- c( "hDC7_mock", "hDC7_ZIKVneg", "hDC7_ZIKVpos",
                         "hDC8_mock", "hDC8_ZIKVneg", "hDC8_ZIKVpos",
                         "hDC10_mock", "hDC10_ZIKVneg", "hDC10_ZIKVpos")
```

``` r
#remove columns with 0 value
my_dat2 <- my_dat2[rowSums(my_dat2) != 0,]
#subtract the 'background' mock from sample
ZIKV_minus_mock <- data.frame(my_dat2[,2:3] - my_dat2$hDC7_mock,
                         my_dat2[,5:6] - my_dat2$hDC8_mock,
                         my_dat2[,8:9] - my_dat2$hDC10_mock)
#ZIKV_minus_mock2 <- ZIKV_minus_mock
```

Add annotation information
==========================

``` r
library("AnnotationDbi")
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    ##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    ##     setdiff, sort, table, tapply, union, unique, unsplit, which,
    ##     which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

``` r
library("org.Hs.eg.db")
```

    ## 

``` r
columns(org.Hs.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
    ##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
    ##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
    ## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
    ## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
    ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    ## [25] "UNIGENE"      "UNIPROT"

``` r
ZIKV_minus_mock$symbol = mapIds(org.Hs.eg.db,
                    keys=rownames(ZIKV_minus_mock), 
                    keytype="REFSEQ",
                    column="SYMBOL",
                    multiVals="first")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
ZIKV_minus_mock$entrez = mapIds(org.Hs.eg.db,
                    keys=rownames(ZIKV_minus_mock),
                    keytype="REFSEQ",
                    column="ENTREZID",
                    multiVals="first")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
ZIKV_minus_mock$name =   mapIds(org.Hs.eg.db,
                    keys=rownames(ZIKV_minus_mock),
                    keytype= "REFSEQ",
                    column= "GENENAME",
                    multiVals="first")
```

    ## 'select()' returned 1:1 mapping between keys and columns

loading libraries for data analysis

``` r
library(pathview)
```

    ## ##############################################################################
    ## Pathview is an open source software package distributed under GNU General
    ## Public License version 3 (GPLv3). Details of GPLv3 is available at
    ## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    ## formally cite the original Pathview paper (not just mention it) in publications
    ## or products. For details, do citation("pathview") within R.
    ## 
    ## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    ## license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ## ##############################################################################

``` r
library(gage)
library(gageData)
```

``` r
data(kegg.sets.hs)
data(sigmet.idx.hs)


kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
```

``` r
eight_neg = ZIKV_minus_mock[,3]
names(eight_neg) = ZIKV_minus_mock$entrez
head(eight_neg)
```

    ##     79630     27125 100125288 105373021     11078 100129387 
    ##     0.125     0.110    -0.586    -0.137     0.434    -0.610

``` r
eight_pos = ZIKV_minus_mock[,4]
names(eight_pos) = ZIKV_minus_mock$entrez
head(eight_pos)
```

    ##     79630     27125 100125288 105373021     11078 100129387 
    ##    -0.367    -0.020     0.348     0.057     0.157    -0.454

``` r
eight_neg_set = gage(eight_neg, gsets= kegg.sets.hs)
eight_pos_set = gage(eight_pos, gsets= kegg.sets.hs)
```

``` r
attributes(eight_neg_set)
```

    ## $names
    ## [1] "greater" "less"    "stats"

``` r
attributes(eight_pos_set)
```

    ## $names
    ## [1] "greater" "less"    "stats"

Display the top results of up and down regulation for both ZIKV+ and ZIKV- samples from hDC8.

``` r
head(eight_neg_set$greater)
```

    ##                                                   p.geomean stat.mean
    ## hsa04623 Cytosolic DNA-sensing pathway         0.0007745945  3.312365
    ## hsa04622 RIG-I-like receptor signaling pathway 0.0014756957  3.053559
    ## hsa04620 Toll-like receptor signaling pathway  0.0022337694  2.894048
    ## hsa04062 Chemokine signaling pathway           0.0029598393  2.773085
    ## hsa04210 Apoptosis                             0.0076929644  2.451983
    ## hsa04380 Osteoclast differentiation            0.0078003812  2.435520
    ##                                                       p.val     q.val
    ## hsa04623 Cytosolic DNA-sensing pathway         0.0007745945 0.1202692
    ## hsa04622 RIG-I-like receptor signaling pathway 0.0014756957 0.1202692
    ## hsa04620 Toll-like receptor signaling pathway  0.0022337694 0.1206135
    ## hsa04062 Chemokine signaling pathway           0.0029598393 0.1206135
    ## hsa04210 Apoptosis                             0.0076929644 0.2119104
    ## hsa04380 Osteoclast differentiation            0.0078003812 0.2119104
    ##                                                set.size         exp1
    ## hsa04623 Cytosolic DNA-sensing pathway               50 0.0007745945
    ## hsa04622 RIG-I-like receptor signaling pathway       66 0.0014756957
    ## hsa04620 Toll-like receptor signaling pathway        96 0.0022337694
    ## hsa04062 Chemokine signaling pathway                187 0.0029598393
    ## hsa04210 Apoptosis                                   84 0.0076929644
    ## hsa04380 Osteoclast differentiation                 126 0.0078003812

``` r
head(eight_pos_set$greater)
```

    ##                                                      p.geomean stat.mean
    ## hsa04623 Cytosolic DNA-sensing pathway            0.0006341159  3.358598
    ## hsa00100 Steroid biosynthesis                     0.0006856734  3.522049
    ## hsa00900 Terpenoid backbone biosynthesis          0.0016587509  3.233243
    ## hsa04620 Toll-like receptor signaling pathway     0.0054810385  2.578627
    ## hsa04610 Complement and coagulation cascades      0.0188261647  2.106488
    ## hsa00260 Glycine, serine and threonine metabolism 0.0193905764  2.114569
    ##                                                          p.val      q.val
    ## hsa04623 Cytosolic DNA-sensing pathway            0.0006341159 0.05588238
    ## hsa00100 Steroid biosynthesis                     0.0006856734 0.05588238
    ## hsa00900 Terpenoid backbone biosynthesis          0.0016587509 0.09012547
    ## hsa04620 Toll-like receptor signaling pathway     0.0054810385 0.22335232
    ## hsa04610 Complement and coagulation cascades      0.0188261647 0.52677733
    ## hsa00260 Glycine, serine and threonine metabolism 0.0193905764 0.52677733
    ##                                                   set.size         exp1
    ## hsa04623 Cytosolic DNA-sensing pathway                  50 0.0006341159
    ## hsa00100 Steroid biosynthesis                           19 0.0006856734
    ## hsa00900 Terpenoid backbone biosynthesis                14 0.0016587509
    ## hsa04620 Toll-like receptor signaling pathway           96 0.0054810385
    ## hsa04610 Complement and coagulation cascades            65 0.0188261647
    ## hsa00260 Glycine, serine and threonine metabolism       30 0.0193905764

``` r
head(eight_neg_set$less)
```

    ##                                                      p.geomean stat.mean
    ## hsa04142 Lysosome                                   0.01936288 -2.078621
    ## hsa00280 Valine, leucine and isoleucine degradation 0.02481877 -1.992037
    ## hsa00640 Propanoate metabolism                      0.02988426 -1.917630
    ## hsa04150 mTOR signaling pathway                     0.04434400 -1.719966
    ## hsa00620 Pyruvate metabolism                        0.04606984 -1.706196
    ## hsa03030 DNA replication                            0.04854634 -1.685488
    ##                                                          p.val     q.val
    ## hsa04142 Lysosome                                   0.01936288 0.6401596
    ## hsa00280 Valine, leucine and isoleucine degradation 0.02481877 0.6401596
    ## hsa00640 Propanoate metabolism                      0.02988426 0.6401596
    ## hsa04150 mTOR signaling pathway                     0.04434400 0.6401596
    ## hsa00620 Pyruvate metabolism                        0.04606984 0.6401596
    ## hsa03030 DNA replication                            0.04854634 0.6401596
    ##                                                     set.size       exp1
    ## hsa04142 Lysosome                                        120 0.01936288
    ## hsa00280 Valine, leucine and isoleucine degradation       43 0.02481877
    ## hsa00640 Propanoate metabolism                            32 0.02988426
    ## hsa04150 mTOR signaling pathway                           51 0.04434400
    ## hsa00620 Pyruvate metabolism                              40 0.04606984
    ## hsa03030 DNA replication                                  36 0.04854634

``` r
head(eight_pos_set$less)
```

    ##                                                p.geomean stat.mean
    ## hsa03010 Ribosome                           1.475930e-12 -7.711653
    ## hsa04810 Regulation of actin cytoskeleton   1.385170e-04 -3.670224
    ## hsa04514 Cell adhesion molecules (CAMs)     4.456362e-03 -2.640795
    ## hsa00190 Oxidative phosphorylation          8.985599e-03 -2.389872
    ## hsa04270 Vascular smooth muscle contraction 9.509604e-03 -2.363213
    ## hsa04144 Endocytosis                        1.143270e-02 -2.284987
    ##                                                    p.val        q.val
    ## hsa03010 Ribosome                           1.475930e-12 2.405766e-10
    ## hsa04810 Regulation of actin cytoskeleton   1.385170e-04 1.128914e-02
    ## hsa04514 Cell adhesion molecules (CAMs)     4.456362e-03 2.421290e-01
    ## hsa00190 Oxidative phosphorylation          8.985599e-03 2.449146e-01
    ## hsa04270 Vascular smooth muscle contraction 9.509604e-03 2.449146e-01
    ## hsa04144 Endocytosis                        1.143270e-02 2.449146e-01
    ##                                             set.size         exp1
    ## hsa03010 Ribosome                                 87 1.475930e-12
    ## hsa04810 Regulation of actin cytoskeleton        200 1.385170e-04
    ## hsa04514 Cell adhesion molecules (CAMs)          125 4.456362e-03
    ## hsa00190 Oxidative phosphorylation               116 8.985599e-03
    ## hsa04270 Vascular smooth muscle contraction      109 9.509604e-03
    ## hsa04144 Endocytosis                             199 1.143270e-02

``` r
pathview(gene.data = eight_pos, pathway.id = "hsa04623", kegg.native = TRUE)
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Info: Working in directory /Users/jvalls/Desktop/BGGN213/ZIKV_DC

    ## Info: Writing image file hsa04623.pathview.png

``` r
pathview(gene.data = eight_neg, pathway.id = "hsa04623", kegg.native = TRUE)
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Info: Working in directory /Users/jvalls/Desktop/BGGN213/ZIKV_DC

    ## Info: Writing image file hsa04623.pathview.png

``` r
# "hsa00061" lipid synt
#  "hsa04979" cholesterol synthesis
```

Combining all three samples into one
====================================

``` r
#intialize with two columns
ZIKV_minus_mock_combine <- data.frame(ZIKV_minus_mock[,1:2])

ZIKV_minus_mock_combine[,1] <- (ZIKV_minus_mock$hDC7_ZIKVneg + ZIKV_minus_mock$hDC8_ZIKVneg + ZIKV_minus_mock$hDC10_ZIKVneg / 3)

ZIKV_minus_mock_combine[,2] <- (ZIKV_minus_mock$hDC7_ZIKVpos + ZIKV_minus_mock$hDC8_ZIKVpos + ZIKV_minus_mock$hDC10_ZIKVpos / 3)

colnames(ZIKV_minus_mock_combine) <- c("neg", "pos")
```

``` r
ZIKV_minus_mock_combine$symbol = mapIds(org.Hs.eg.db,
                    keys=rownames(ZIKV_minus_mock_combine), 
                    keytype="REFSEQ",
                    column="SYMBOL",
                    multiVals="first")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
ZIKV_minus_mock_combine$entrez = mapIds(org.Hs.eg.db,
                    keys=rownames(ZIKV_minus_mock_combine),
                    keytype="REFSEQ",
                    column="ENTREZID",
                    multiVals="first")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
ZIKV_minus_mock_combine$name =   mapIds(org.Hs.eg.db,
                    keys=rownames(ZIKV_minus_mock_combine),
                    keytype= "REFSEQ",
                    column= "GENENAME",
                    multiVals="first")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
combine_neg = ZIKV_minus_mock_combine[,1]
names(combine_neg) = ZIKV_minus_mock_combine$entrez
head(combine_neg)
```

    ##      79630      27125  100125288  105373021      11078  100129387 
    ##  0.4690000  0.3036667 -1.3776667 -0.6780000  0.7503333 -0.5343333

``` r
combine_pos = ZIKV_minus_mock_combine[,2]
names(combine_pos) = ZIKV_minus_mock_combine$entrez
head(combine_pos)
```

    ##       79630       27125   100125288   105373021       11078   100129387 
    ## -0.71866667  0.16366667  0.52000000 -0.14400000  0.30200000  0.05566667

``` r
combine_neg_set = gage(combine_neg, gsets= kegg.sets.hs)
combine_pos_set = gage(combine_pos, gsets= kegg.sets.hs)
```

``` r
head(combine_neg_set$greater)
```

    ##                                                   p.geomean stat.mean
    ## hsa04062 Chemokine signaling pathway           2.215788e-05  4.149554
    ## hsa04620 Toll-like receptor signaling pathway  1.447256e-04  3.726229
    ## hsa04630 Jak-STAT signaling pathway            1.450151e-04  3.672278
    ## hsa04623 Cytosolic DNA-sensing pathway         1.632705e-04  3.806845
    ## hsa04622 RIG-I-like receptor signaling pathway 1.686925e-04  3.721775
    ## hsa03050 Proteasome                            6.926793e-04  3.319576
    ##                                                       p.val       q.val
    ## hsa04062 Chemokine signaling pathway           2.215788e-05 0.003611734
    ## hsa04620 Toll-like receptor signaling pathway  1.447256e-04 0.005499375
    ## hsa04630 Jak-STAT signaling pathway            1.450151e-04 0.005499375
    ## hsa04623 Cytosolic DNA-sensing pathway         1.632705e-04 0.005499375
    ## hsa04622 RIG-I-like receptor signaling pathway 1.686925e-04 0.005499375
    ## hsa03050 Proteasome                            6.926793e-04 0.018817789
    ##                                                set.size         exp1
    ## hsa04062 Chemokine signaling pathway                187 2.215788e-05
    ## hsa04620 Toll-like receptor signaling pathway        96 1.447256e-04
    ## hsa04630 Jak-STAT signaling pathway                 136 1.450151e-04
    ## hsa04623 Cytosolic DNA-sensing pathway               50 1.632705e-04
    ## hsa04622 RIG-I-like receptor signaling pathway       66 1.686925e-04
    ## hsa03050 Proteasome                                  42 6.926793e-04

``` r
head(combine_pos_set$greater)
```

    ##                                                   p.geomean stat.mean
    ## hsa00100 Steroid biosynthesis                  0.0004424800  3.742330
    ## hsa04623 Cytosolic DNA-sensing pathway         0.0006332482  3.365687
    ## hsa00900 Terpenoid backbone biosynthesis       0.0019228839  3.199465
    ## hsa04620 Toll-like receptor signaling pathway  0.0021891417  2.899290
    ## hsa04622 RIG-I-like receptor signaling pathway 0.0102070424  2.354426
    ## hsa04610 Complement and coagulation cascades   0.0215398534  2.050933
    ##                                                       p.val      q.val
    ## hsa00100 Steroid biosynthesis                  0.0004424800 0.05160972
    ## hsa04623 Cytosolic DNA-sensing pathway         0.0006332482 0.05160972
    ## hsa00900 Terpenoid backbone biosynthesis       0.0019228839 0.08920752
    ## hsa04620 Toll-like receptor signaling pathway  0.0021891417 0.08920752
    ## hsa04622 RIG-I-like receptor signaling pathway 0.0102070424 0.33274958
    ## hsa04610 Complement and coagulation cascades   0.0215398534 0.58516602
    ##                                                set.size         exp1
    ## hsa00100 Steroid biosynthesis                        19 0.0004424800
    ## hsa04623 Cytosolic DNA-sensing pathway               50 0.0006332482
    ## hsa00900 Terpenoid backbone biosynthesis             14 0.0019228839
    ## hsa04620 Toll-like receptor signaling pathway        96 0.0021891417
    ## hsa04622 RIG-I-like receptor signaling pathway       66 0.0102070424
    ## hsa04610 Complement and coagulation cascades         65 0.0215398534

``` r
head(combine_neg_set$less)
```

    ##                                                      p.geomean stat.mean
    ## hsa00450 Selenocompound metabolism                  0.02688177 -2.002763
    ## hsa03410 Base excision repair                       0.03801122 -1.809471
    ## hsa04150 mTOR signaling pathway                     0.04448284 -1.717840
    ## hsa03430 Mismatch repair                            0.04573443 -1.735950
    ## hsa00280 Valine, leucine and isoleucine degradation 0.04712889 -1.693295
    ## hsa00640 Propanoate metabolism                      0.04881523 -1.682239
    ##                                                          p.val     q.val
    ## hsa00450 Selenocompound metabolism                  0.02688177 0.8765422
    ## hsa03410 Base excision repair                       0.03801122 0.8765422
    ## hsa04150 mTOR signaling pathway                     0.04448284 0.8765422
    ## hsa03430 Mismatch repair                            0.04573443 0.8765422
    ## hsa00280 Valine, leucine and isoleucine degradation 0.04712889 0.8765422
    ## hsa00640 Propanoate metabolism                      0.04881523 0.8765422
    ##                                                     set.size       exp1
    ## hsa00450 Selenocompound metabolism                        17 0.02688177
    ## hsa03410 Base excision repair                             33 0.03801122
    ## hsa04150 mTOR signaling pathway                           51 0.04448284
    ## hsa03430 Mismatch repair                                  23 0.04573443
    ## hsa00280 Valine, leucine and isoleucine degradation       43 0.04712889
    ## hsa00640 Propanoate metabolism                            32 0.04881523

``` r
head(combine_pos_set$less)
```

    ##                                                     p.geomean stat.mean
    ## hsa03010 Ribosome                                3.817331e-10 -6.712332
    ## hsa04810 Regulation of actin cytoskeleton        4.099960e-04 -3.372711
    ## hsa04142 Lysosome                                9.515966e-03 -2.362643
    ## hsa04914 Progesterone-mediated oocyte maturation 1.276184e-02 -2.253707
    ## hsa04070 Phosphatidylinositol signaling system   1.276496e-02 -2.257045
    ## hsa03040 Spliceosome                             1.528565e-02 -2.177862
    ##                                                         p.val        q.val
    ## hsa03010 Ribosome                                3.817331e-10 6.222249e-08
    ## hsa04810 Regulation of actin cytoskeleton        4.099960e-04 3.341468e-02
    ## hsa04142 Lysosome                                9.515966e-03 3.653329e-01
    ## hsa04914 Progesterone-mediated oocyte maturation 1.276184e-02 3.653329e-01
    ## hsa04070 Phosphatidylinositol signaling system   1.276496e-02 3.653329e-01
    ## hsa03040 Spliceosome                             1.528565e-02 3.653329e-01
    ##                                                  set.size         exp1
    ## hsa03010 Ribosome                                      87 3.817331e-10
    ## hsa04810 Regulation of actin cytoskeleton             200 4.099960e-04
    ## hsa04142 Lysosome                                     120 9.515966e-03
    ## hsa04914 Progesterone-mediated oocyte maturation       85 1.276184e-02
    ## hsa04070 Phosphatidylinositol signaling system         74 1.276496e-02
    ## hsa03040 Spliceosome                                  127 1.528565e-02
