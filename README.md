## Dependancies and setup

Set up colors and themes.

    theme_set(theme_bw())
    theme_update(plot.title = element_text(hjust = 0.5))
    my_colors = colorRampPalette(c("violet", "white", "darkorange"))(n=599)
    my_colors_2 = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(n=29)

## Load the bladder cancer data

Load the expression data from 57 bladder samples. Samples were processed
in 5 batches. [Original
study](https://cancerres.aacrjournals.org/content/canres/64/11/4040.full.pdf)

    data(bladderdata)
    phenotypes <- pData(bladderEset)
    expression_data <- exprs(bladderEset)
    phenotypes

    ## Warning in knit_print.huxtable(ht): Unrecognized output format "markdown". Using `to_screen` to print huxtables.
    ## Set options("huxtable.knitr_output_format") manually to "latex", "html", "rtf", "docx", "pptx", "md" or "screen".

                     ┌────────────────────────────────────┐
                     │ sample   outcome    batch   cancer │
                     ├────────────────────────────────────┤
                     │      1   Normal         3   Normal │
                     │      2   Normal         2   Normal │
                     │      3   Normal         2   Normal │
                     │      4   Normal         3   Normal │
                     │      5   Normal         3   Normal │
                     │      6   Normal         3   Normal │
                     │      7   Normal         2   Normal │
                     │      8   Normal         2   Normal │
                     │      9   sTCC+CIS       5   Cancer │
                     │     10   sTCC-CIS       2   Cancer │
                     │     11   sTCC-CIS       5   Cancer │
                     │     12   sTCC-CIS       2   Cancer │
                     │     13   sTCC+CIS       5   Cancer │
                     │     14   sTCC-CIS       2   Cancer │
                     │     15   sTCC+CIS       5   Cancer │
                     │     16   sTCC+CIS       5   Cancer │
                     │     17   sTCC-CIS       2   Cancer │
                     │     18   mTCC           1   Cancer │
                     │     19   sTCC+CIS       5   Cancer │
                     │     20   mTCC           1   Cancer │
                     │     21   mTCC           2   Cancer │
                     │     22   mTCC           1   Cancer │
                     │     23   sTCC-CIS       2   Cancer │
                     │     24   sTCC+CIS       5   Cancer │
                     │     25   sTCC-CIS       2   Cancer │
                     │     26   sTCC-CIS       2   Cancer │
                     │     27   sTCC+CIS       5   Cancer │
                     │     28   mTCC           1   Cancer │
                     │     29   mTCC           1   Cancer │
                     │     30   sTCC-CIS       2   Cancer │
                     │     31   mTCC           1   Cancer │
                     │     32   mTCC           1   Cancer │
                     │     33   mTCC           1   Cancer │
                     │     34   sTCC+CIS       5   Cancer │
                     │     35   mTCC           1   Cancer │
                     │     36   sTCC-CIS       2   Cancer │
                     │     37   sTCC-CIS       2   Cancer │
                     │     38   sTCC-CIS       2   Cancer │
                     │     39   sTCC-CIS       2   Cancer │
                     │     40   mTCC           1   Cancer │
                     │     41   sTCC+CIS       5   Cancer │
                     │     42   sTCC+CIS       5   Cancer │
                     │     43   sTCC+CIS       5   Cancer │
                     │     44   sTCC-CIS       2   Cancer │
                     │     45   sTCC-CIS       5   Cancer │
                     │     46   mTCC           1   Cancer │
                     │     47   sTCC-CIS       5   Cancer │
                     │     48   sTCC+CIS       5   Cancer │
                     │     49   Biopsy         4   Biopsy │
                     │     50   Biopsy         4   Biopsy │
                     │     51   Biopsy         5   Biopsy │
                     │     52   Biopsy         5   Biopsy │
                     │     53   Biopsy         5   Biopsy │
                     │     54   Biopsy         5   Biopsy │
                     │     55   Biopsy         4   Biopsy │
                     │     56   Biopsy         4   Biopsy │
                     │     57   Biopsy         4   Biopsy │
                     └────────────────────────────────────┘

Column names: sample, outcome, batch, cancer

    class(expression_data)

    ## [1] "matrix" "array"

Look for N/A nad NaN values in the datasets

    print("Expression data:")

    ## [1] "Expression data:"

    sum(is.na(expression_data))

    ## [1] 0

    print("Metadata:")

    ## [1] "Metadata:"

    sum(is.na(phenotypes))

    ## [1] 0

Plot mean vs variance

    row.variances <- apply(expression_data, 1, function(x) var(x))
    row.means <- apply(expression_data, 1, function(x) mean(x))
    ggplot(data.frame(row.means, row.variances), aes(x=row.means, y = row.variances)) + geom_point() + geom_smooth(method="lm", col = "lightcoral") + labs(title = "Mean vs. Variance relationship") 

    ## `geom_smooth()` using formula = 'y ~ x'

![](homework_files/figure-markdown_strict/plot%20the%20mean%20vs.%20variance-1.png)

Normalize the data

    expression_data <- expression_data[row.variances < 5, ]
    expression_data <- log2(expression_data)
    row.variances <- apply(expression_data, 1, function(x) var(x))
    row.means <- apply(expression_data, 1, function(x) mean(x))
    ggplot(data.frame(row.means, row.variances), aes(x=row.means, y = row.variances)) + geom_point() + geom_smooth(method="lm", col = "lightcoral") + labs(title = "Mean vs. Variance relationship") 

    ## `geom_smooth()` using formula = 'y ~ x'

![](homework_files/figure-markdown_strict/normalize%20the%20data%20and%20plot%20mean%20vs.%20variance-1.png)

## Homework problem 1

Create the table

    #batch.table <- unique(phenotypes$batch) %>% map_dfc(setNames, object = phenotypes$outcome)
    batch.table <- table(phenotypes$cancer, phenotypes$batch, dnn = c("Batch", "Outcome"))
    batch.table <- as_huxtable(batch.table)
    caption(batch.table) <- "Number of certain outcome by batch number"
    bottom_border(batch.table)[1, 1:ncol(batch.table)] <- 0.5
    batch.table <- set_background_color(batch.table, evens, everywhere, "gray95")
    batch.table <- set_bold(batch.table, row = 1, col = everywhere)
    batch.table <- set_bold(batch.table, col = 1, row = everywhere)
    #quick_pdf(batch.table, file = "Guzek_problem1.pdf") 
    print_screen(batch.table, file = "Guzek_problem1.pdf") 

    ##                    Number of certain outcome by batch number                    
    ##                                  1     2     3     4   5    
    ##                     ────────────────────────────────────────
    ##                       Biopsy     0     0     0     5   4    
    ##                       Cancer    11    14     0     0   15   
    ##                       Normal     0     4     4     0   0    
    ## 
    ## Column names: rownames, 1, 2, 3, 4, 5

## Homework problem 2

Sort the data by batches

    # Clumsy way to sort somples by batches from which they come from
    ettiquetted_edf <- arrange(data.frame(phenotypes, t(expression_data)), batch)
    expression_data <- t(as.matrix(arrange(data.frame(t(expression_data)), phenotypes$batch)))
    expression_data[1,]

    ## GSM71037.CEL GSM71039.CEL GSM71041.CEL GSM71047.CEL GSM71048.CEL GSM71050.CEL 
    ##     3.329732     3.221615     3.363534     3.227783     3.369852     3.344546 
    ## GSM71051.CEL GSM71052.CEL GSM71054.CEL GSM71060.CEL GSM71066.CEL GSM71020.CEL 
    ##     3.185236     3.215314     3.309036     3.340781     3.363495     3.109034 
    ## GSM71021.CEL GSM71025.CEL GSM71026.CEL GSM71029.CEL GSM71031.CEL GSM71033.CEL 
    ##     3.134095     3.187140     3.126909     3.346049     3.349686     3.266322 
    ## GSM71036.CEL GSM71040.CEL GSM71042.CEL GSM71044.CEL GSM71045.CEL GSM71049.CEL 
    ##     3.411365     3.349954     3.367420     3.267263     3.391411     3.293280 
    ## GSM71055.CEL GSM71056.CEL GSM71058.CEL GSM71059.CEL GSM71064.CEL GSM71019.CEL 
    ##     3.393006     3.380965     3.309156     3.398594     3.321140     3.338449 
    ## GSM71022.CEL GSM71023.CEL GSM71024.CEL GSM71069.CEL GSM71070.CEL GSM71075.CEL 
    ##     3.209230     3.358515     3.325262     3.323230     3.309300     3.248161 
    ## GSM71076.CEL GSM71077.CEL GSM71028.CEL GSM71030.CEL GSM71032.CEL GSM71034.CEL 
    ##     3.254105     3.176186     3.293263     3.073975     3.348261     3.368876 
    ## GSM71035.CEL GSM71038.CEL GSM71043.CEL GSM71046.CEL GSM71053.CEL GSM71061.CEL 
    ##     3.167771     3.288066     3.220811     3.438364     3.268539     3.279898 
    ## GSM71062.CEL GSM71063.CEL GSM71065.CEL GSM71067.CEL GSM71068.CEL GSM71071.CEL 
    ##     3.378778     3.131505     3.291425     3.409690     3.403662     3.185360 
    ## GSM71072.CEL GSM71073.CEL GSM71074.CEL 
    ##     3.171836     3.069713     3.154367

Prepare the cleaned dataset using ComBat

    batch <- phenotypes$batch
    combat_expression_data <- ComBat(dat = expression_data,
                                     batch = batch, 
                                     mod = model.matrix(~1, data = phenotypes),
                                     par.prior = TRUE,
                                     prior.plots = TRUE)

    ## Found5batches

    ## Adjusting for0covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding parametric adjustments

    ## Adjusting the Data

![](homework_files/figure-markdown_strict/ComBat%20on%20bladderbatch-1.png)
Take a look at data cleaned by ComBat

    head(combat_expression_data[,1:6])

    ##            GSM71037.CEL GSM71039.CEL GSM71041.CEL GSM71047.CEL GSM71048.CEL
    ## X1007_s_at     3.321690     3.221841     3.361954     3.183786     3.375958
    ## X1053_at       2.440131     2.372637     2.422623     2.520566     2.370984
    ## X117_at        2.642628     2.655768     2.595562     2.708171     2.509199
    ## X121_at        3.081410     2.990323     3.025190     3.175983     3.172938
    ## X1255_g_at     1.994907     2.023182     1.959205     2.058130     1.937482
    ## X1294_at       2.844929     2.935953     2.945138     2.865110     2.993894
    ##            GSM71050.CEL
    ## X1007_s_at     3.341727
    ## X1053_at       2.529036
    ## X117_at        2.584665
    ## X121_at        3.005630
    ## X1255_g_at     1.950303
    ## X1294_at       3.075656

Create the heatmaps

    #png("Guzek_problem2a.png", width = 900, height = 900)
    heatmap.2(as.matrix(expression_data),
              main = "Bladder Cancer Data Sorted by Bath, Before ComBat",
              notecol = "black",
              Colv = FALSE,
              density.info = "none",
              trace = "none",
              margins = c(12,9),
              col = my_colors,
              dendrogram = "none",
              scale = "row",
              ColSideColors = palette("Set2")[as.factor(ettiquetted_edf$batch)]
    )

![](homework_files/figure-markdown_strict/heatmap%20before%20ComBat-1.png)

    #dev.off()

    #png("Guzek_problem2b.png", width = 900, height = 900)
    heatmap.2(as.matrix(combat_expression_data),
              main = "Bladder Cancer Data Sorted by Batch, After Combat",
              notecol = "black",
              Colv = FALSE,
              density.info = "none",
              trace = "none",
              margins = c(12,9),
              col = my_colors,
              dendrogram = "none",
              scale = "row",
              ColSideColors = palette("Set2")[as.factor(ettiquetted_edf$batch)]
    )

![](homework_files/figure-markdown_strict/heatmap%20after%20ComBat-1.png)

    #dev.off()

# Homework problem 3

    correlation.matrix <- apply(expression_data, 2, cor, expression_data)
    rownames(correlation.matrix) <- colnames(correlation.matrix)
    #png("Guzek_problem3a.png", width = 1000, height = 1000)
    heatmap.2(correlation.matrix,
              main = "Pairwise Pearson Comparison of Expression Data, Samples Grouped by Batch",
              notecol = "black",
              density.info = "none",
              Colv = F,
              Rowv = F,
              trace = "none",
              colsep = 1:ncol(correlation.matrix),
              rowsep = 1:ncol(correlation.matrix),
              sepcolor = "gray",
              margins = c(9, 6),
              dendrogram = "none",
              col = my_colors_2,
              RowSideColors = palette("Set2")[as.factor(ettiquetted_edf$batch)],
              ColSideColors =  palette("Set2")[as.factor(ettiquetted_edf$batch)]
    )

![](homework_files/figure-markdown_strict/Correlation%20heatmaps-1.png)

    #dev.off()

    ettiquetted_edf = arrange(ettiquetted_edf, outcome)
    expression_data = t(as.matrix(arrange(data.frame(t(expression_data)), phenotypes$outcome)))

    correlation.matrix <- apply(expression_data, 2, cor, expression_data)
    rownames(correlation.matrix) <- colnames(correlation.matrix)
    #png("Guzek_problem3b.png", width = 1000, height = 1000)
    heatmap.2(correlation.matrix,
              main = "Pairwise Pearson Comparison of Expression Data, Samples Grouped by Outcome",
              notecol = "black",
              density.info = "none",
              Colv = F,
              Rowv = F,
              trace = "none",
              colsep = 1:ncol(correlation.matrix),
              rowsep = 1:ncol(correlation.matrix),
              sepcolor = "gray",
              margins = c(9, 6),
              dendrogram = "none",
              col = my_colors_2,
              RowSideColors = palette("Set2")[as.factor(ettiquetted_edf$outcome)],
              ColSideColors =  palette("Set2")[as.factor(ettiquetted_edf$outcome)]
    )

![](homework_files/figure-markdown_strict/Correlation%20heatmaps-2.png)

    #dev.off()

## Homework problem 4

Import the bottomly datasest

    if (file_test("-f", "./bottomly.Rdata")) {
      load(file = "./bottomly.Rdata")
    } else {
      con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
      load(file = con)
      close(con)
      save(bottomly.eset, file = "bottomly.Rdata")
    }

Get the expression data

    bottomly_edata <- as.matrix(exprs(bottomly.eset))
    bottomly_edata <- bottomly_edata[rowMeans(bottomly_edata) > 10, ]
    bottomly_edata <- log2(as.matrix(bottomly_edata) + 1)
    bottomly_metadata <- pData(bottomly.eset)
    dim(bottomly_edata)

    ## [1] 8544   21

    bottomly_edata[1:10,1:5]

    ##                    SRX033480 SRX033488 SRX033481 SRX033489 SRX033482
    ## ENSMUSG00000000001  8.531381  9.541097  8.169925  9.588715  8.447083
    ## ENSMUSG00000000056  4.459432  5.554589  4.392317  5.209453  3.700440
    ## ENSMUSG00000000058  4.000000  5.459432  3.700440  5.129283  3.906891
    ## ENSMUSG00000000078  9.016808  9.773139  8.413628  9.668885  8.566054
    ## ENSMUSG00000000088  8.098032  9.611025  8.108524  9.543032  8.238405
    ## ENSMUSG00000000093  3.584963  4.643856  3.807355  5.129283  4.807355
    ## ENSMUSG00000000120  4.321928  4.906891  4.807355  6.794416  4.807355
    ## ENSMUSG00000000125  4.321928  4.643856  4.321928  5.426265  4.523562
    ## ENSMUSG00000000126  6.614710  7.607330  5.700440  6.672425  6.087463
    ## ENSMUSG00000000127  6.768184  7.467606  6.392317  7.383704  6.303781

    bottomly_metadata

    ## Warning in knit_print.huxtable(ht): Unrecognized output format "markdown". Using `to_screen` to print huxtables.
    ## Set options("huxtable.knitr_output_format") manually to "latex", "html", "rtf", "docx", "pptx", "md" or "screen".

      ┌──────────────────────────────────────────────────────────────────┐
      │ sample.id   num.tech.rep   strain     experiment.n   lane.number │
      │                        s                     umber               │
      ├──────────────────────────────────────────────────────────────────┤
      │ SRX033480              1   C57BL/6J              6             1 │
      │ SRX033488              1   C57BL/6J              7             1 │
      │ SRX033481              1   C57BL/6J              6             2 │
      │ SRX033489              1   C57BL/6J              7             2 │
      │ SRX033482              1   C57BL/6J              6             3 │
      │ SRX033490              1   C57BL/6J              7             3 │
      │ SRX033483              1   C57BL/6J              6             5 │
      │ SRX033476              1   C57BL/6J              4             6 │
      │ SRX033478              1   C57BL/6J              4             7 │
      │ SRX033479              1   C57BL/6J              4             8 │
      │ SRX033472              1   DBA/2J                4             1 │
      │ SRX033473              1   DBA/2J                4             2 │
      │ SRX033474              1   DBA/2J                4             3 │
      │ SRX033475              1   DBA/2J                4             5 │
      │ SRX033491              1   DBA/2J                7             5 │
      │ SRX033484              1   DBA/2J                6             6 │
      │ SRX033492              1   DBA/2J                7             6 │
      │ SRX033485              1   DBA/2J                6             7 │
      │ SRX033493              1   DBA/2J                7             7 │
      │ SRX033486              1   DBA/2J                6             8 │
      │ SRX033494              1   DBA/2J                7             8 │
      └──────────────────────────────────────────────────────────────────┘

Column names: sample.id, num.tech.reps, strain, experiment.number,
lane.number

    mod = lm(t(bottomly_edata) ~ as.factor(bottomly_metadata$strain) + as.factor(bottomly_metadata$experiment.number))
    tidy_mod <- tidy(mod)
    head(tidy_mod)

    ## Warning in knit_print.huxtable(ht): Unrecognized output format "markdown". Using `to_screen` to print huxtables.
    ## Set options("huxtable.knitr_output_format") manually to "latex", "html", "rtf", "docx", "pptx", "md" or "screen".

┌───────────────────────────────────────────────────────────────────────┐
│ response term estimate std.error statistic p.value │
├───────────────────────────────────────────────────────────────────────┤
│ ENSMUSG000 (Intercept 8.63   0.135 64     1.05e-21 │ │ 00000001 ) │ │
ENSMUSG000 as.factor( -0.0207 0.131 -0.158 0.876    │ │ 00000001
bottomly\_m │ │
etadata*s**t*││*r**a**i**n*)*D**B**A*/2││*J*││*E**N**S**M**U**S**G*000*a**s*.*f**a**c**t**o**r*( − 0.0560.16 − 0.3510.73││00000001*b**o**t**t**o**m**l**y*<sub>*m*</sub>││*e**t**a**d**a**t**a*ex
│ │ periment.n │ │ umber)6 │ │ ENSMUSG000 as.factor( 1.01   0.159 6.36 
7.1e-06  │ │ 00000001 bottomly\_m │ │
etadata*e**x*││*p**e**r**i**m**e**n**t*.*n*││*u**m**b**e**r*)7││*E**N**S**M**U**S**G*000(*I**n**t**e**r**c**e**p**t*5.20.14535.91.78*e*−17││00000056)││*E**N**S**M**U**S**G*000*a**s*.*f**a**c**t**o**r*(0.04360.1410.310.76││00000056*b**o**t**t**o**m**l**y*<sub>*m*</sub>││*e**t**a**d**a**t**a*st
│ │ rain)DBA/2 │ │ J │
└───────────────────────────────────────────────────────────────────────┘

Column names: response, term, estimate, std.error, statistic, p.value

    names(tidy_mod)

    ## [1] "response"  "term"      "estimate"  "std.error" "statistic" "p.value"

    png("Guzek_problem4a.png", width = 800, height = 700)
    coeff = ggplot(tidy_mod %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J")) + geom_histogram(aes(x=estimate), bins = 100, fill = "darkorange")
    p_val = ggplot(tidy_mod %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")
    figure <- ggarrange(coeff, p_val, labels = c("Histogram of coefficients of uncorrected data", "Histogram of p-values of uncorrected data"), ncol = 1, nrow = 2) + labs(title = "Guzek_problem4a")
    annotate_figure(figure, top = text_grob("Guzek_problem4a"))
    dev.off()

    ## png 
    ##   2

    figure

![](homework_files/figure-markdown_strict/Bottomly%20lm%20histogram-1.png)

    bottomly_batch <- bottomly_metadata$experiment.number
    combat_bottomly_edata <- ComBat(dat = bottomly_edata,
                                    batch = bottomly_metadata$experiment.number,
                                    mod = model.matrix(~1, data = bottomly_metadata),
                                    par.prior = TRUE,
                                    prior.plots = TRUE
    )

    ## Found3batches

    ## Adjusting for0covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding parametric adjustments

    ## Adjusting the Data

![](homework_files/figure-markdown_strict/ComBat%20on%20bottomly%20data-1.png)

Fit linear model on combat corrected bottomly data.

    mod_combat <- lm(t(combat_bottomly_edata) ~ as.factor(bottomly_metadata$strain))
    tidy_mod_combat <- tidy(mod_combat)
    head(tidy_mod_combat)

    ## Warning in knit_print.huxtable(ht): Unrecognized output format "markdown". Using `to_screen` to print huxtables.
    ## Set options("huxtable.knitr_output_format") manually to "latex", "html", "rtf", "docx", "pptx", "md" or "screen".

┌───────────────────────────────────────────────────────────────────────┐
│ response term estimate std.error statistic p.value │
├───────────────────────────────────────────────────────────────────────┤
│ ENSMUSG000 (Intercept 8.93    0.0762 117      1.24e-28 │ │ 00000001 )
│ │ ENSMUSG000 as.factor( -0.00286 0.105  -0.0272 0.979    │ │ 00000001
bottomly\_m │ │
etadata*s**t*││*r**a**i**n*)*D**B**A*/2││*J*││*E**N**S**M**U**S**G*000(*I**n**t**e**r**c**e**p**t*5.070.090755.91.52*e*−22││00000056)││*E**N**S**M**U**S**G*000*a**s*.*f**a**c**t**o**r*(0.03130.1250.250.806││00000056*b**o**t**t**o**m**l**y*<sub>*m*</sub>││*e**t**a**d**a**t**a*st
│ │ rain)DBA/2 │ │ J │ │ ENSMUSG000 (Intercept 4.65    0.103  45.2   
8.46e-21 │ │ 00000058 ) │ │ ENSMUSG000 as.factor( -0.119   0.142 
-0.84   0.411    │ │ 00000058 bottomly\_m │ │ etadata$st │ │ rain)DBA/2
│ │ J │
└───────────────────────────────────────────────────────────────────────┘

Column names: response, term, estimate, std.error, statistic, p.value

    names(tidy_mod_combat)

    ## [1] "response"  "term"      "estimate"  "std.error" "statistic" "p.value"

    png("Guzek_problem4b.png", width = 800, height = 700)
    coeff = ggplot(tidy_mod_combat %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J")) + geom_histogram(aes(x=estimate), bins = 100, fill = "darkorange")
    p_val = ggplot(tidy_mod_combat %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")
    figure <- ggarrange(coeff, p_val, labels = c("Histogram of coeffecients of corrected data", "Histogram of p-values of corrected data"), ncol = 1, nrow = 2) + labs(title = "Guzek_problem4b")
    annotate_figure(figure, top = text_grob("Guzek_problem4b"))
    dev.off()

    ## png 
    ##   2

    figure

![](homework_files/figure-markdown_strict/histogram%20of%20lm%20coefficients%20and%20pvalues-1.png)

Compare the coefficients of uncorrected and corrected data

    est_compare <- tibble(
      Uncorrected = tidy_mod %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("estimate") %>% unlist,
      Corrected = tidy_mod_combat %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("estimate") %>% unlist)

    png("Guzek_problem4c.png", width = 800, height = 500)
    figure <- ggplot(est_compare, aes(x = Uncorrected, y = Corrected)) +
         geom_point(col = "darkgrey", alpha = .5, size = .5) + geom_abline(intercept = 0, slope = 1, col = "darkorange") + geom_smooth(method = "lm") + labs(title = "Guzek_problem4c")
    figure

    ## `geom_smooth()` using formula = 'y ~ x'

    dev.off()

    ## png 
    ##   2

    figure

    ## `geom_smooth()` using formula = 'y ~ x'

![](homework_files/figure-markdown_strict/compare%20coefficients%20of%20lm%20and%20ComBat-1.png)

## Homework problem 5

Find a dimension of surrogate variables

    mod = model.matrix(~as.factor(strain), data = bottomly_metadata)
    mod0 = model.matrix(~1, data = bottomly_metadata)
    num.sv(bottomly_edata, mod, method = "be")

    ## [1] 2

    num.sv(bottomly_edata, mod, method = "leek")

    ## [1] 2

    sva.output <- sva(bottomly_edata, mod, mod0, n.sv=num.sv(bottomly_edata, mod, method = "leek"))

    ## Number of significant surrogate variables is:  2 
    ## Iteration (out of 5 ):1  2  3  4  5

Take a look at the data

    head(sva.output$sv)

    ##            [,1]       [,2]
    ## [1,]  0.2726372 0.26961200
    ## [2,] -0.2092806 0.11715634
    ## [3,]  0.3344130 0.29094759
    ## [4,] -0.2473647 0.10943716
    ## [5,]  0.2788337 0.19859530
    ## [6,] -0.2951497 0.07796639

Look at the relation between the surrogate variables and experiment
number

    summary(lm(sva.output$sv ~ bottomly_metadata$experiment.number))

    ## Response Y1 :
    ## 
    ## Call:
    ## lm(formula = Y1 ~ bottomly_metadata$experiment.number)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.22961 -0.15941 -0.05361  0.16796  0.36437 
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                          0.50935    0.19982   2.549   0.0196 *
    ## bottomly_metadata$experiment.number -0.08988    0.03444  -2.610   0.0172 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1968 on 19 degrees of freedom
    ## Multiple R-squared:  0.2639, Adjusted R-squared:  0.2252 
    ## F-statistic: 6.812 on 1 and 19 DF,  p-value: 0.01721
    ## 
    ## 
    ## Response Y2 :
    ## 
    ## Call:
    ## lm(formula = Y2 ~ bottomly_metadata$experiment.number)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.1409 -0.1167 -0.0700  0.1540  0.2463 
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                         -0.75852    0.14998  -5.058 6.99e-05 ***
    ## bottomly_metadata$experiment.number  0.13386    0.02585   5.179 5.34e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1477 on 19 degrees of freedom
    ## Multiple R-squared:  0.5853, Adjusted R-squared:  0.5635 
    ## F-statistic: 26.82 on 1 and 19 DF,  p-value: 5.344e-05

Visualize the SVs

    sva.batch <- tibble(SV1 = sva.output$sv[,1],
                      SV2 = sva.output$sv[,2],
                      batch = as.factor(bottomly_metadata$experiment.number),
                      strain = as.factor(bottomly_metadata$strain),
                      lane = as.factor(bottomly_metadata$lane.number))

    batch.plot <- ggplot(sva.batch) + geom_point(aes(x = SV1, y = SV2, col = batch))
    batch.plot

![](homework_files/figure-markdown_strict/sva%20vis%201-1.png)

    strain.plot <- ggplot(sva.batch) + geom_point(aes(x = SV1, y = SV2, col = strain))
    strain.plot

![](homework_files/figure-markdown_strict/sva%20vis%202-1.png)

    lane.plot <- ggplot(sva.batch) + geom_point(aes(x = SV1, y = SV2, col = lane))
    lane.plot

![](homework_files/figure-markdown_strict/sva%20vis%203-1.png)

    sva.batch <- tibble(SV1 = sva.output$sv[,1],
                        SV2 = sva.output$sv[,2],
                        batch = as.factor(bottomly_metadata$experiment.number))
    sva.batch.gather <- gather(sva.batch, "sv", "value", -batch)

    ggplot(sva.batch.gather) + geom_violin(aes(x = batch, y = value)) + facet_wrap(~sv, ncol = 1) + geom_jitter(aes(x = batch, y = value, col = batch))

![](homework_files/figure-markdown_strict/sva%20vis%20violin-1.png)

    # Add the surrogate variables to the model matrix
    mod_sva = lm(t(bottomly_edata) ~ as.factor(bottomly_metadata$strain) + sva.output$sv)
    tidy_mod_sva <- tidy(mod_sva)

    est_compare <- tibble(
      LinearModel = tidy_mod %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("estimate") %>% unlist,
      ComBat = tidy_mod_combat %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("estimate") %>% unlist,
      SVA = tidy_mod_sva %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("estimate") %>% unlist)

    png("Guzek_problem5a.png", width = 1500, height = 500)
    lmcb <- ggplot(est_compare, aes(x=LinearModel, y=ComBat)) +
         geom_point(col = "darkgrey", alpha = .5, size = .5) + geom_abline(intercept = 0, slope = 1, col="darkorange") + geom_smooth(method = "lm")

    lmsva <- ggplot(est_compare, aes(x=LinearModel, y=SVA)) +
         geom_point(col = "darkgrey", alpha = .5, size = .5) + geom_abline(intercept = 0, slope = 1, col="darkorange") + geom_smooth(method = "lm")

    cbsva <- ggplot(est_compare, aes(x=ComBat, y=SVA)) +
         geom_point(col = "darkgrey", alpha = .5, size = .5) + geom_abline(intercept = 0, slope = 1, col="darkorange") + geom_smooth(method = "lm")
    figure <- ggarrange(lmcb, lmsva, cbsva, labels = c("Linear Model vs. ComBat", "Linear Model vs. SVA", "ComBat vs. SVA"), nrow = 1, ncol = 3) + labs(title = "Guzek_problem5a")

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    annotate_figure(figure, top = text_grob("Guzek_problem5a"))
    dev.off()

    ## png 
    ##   2

    figure

![](homework_files/figure-markdown_strict/visualize%20the%20comparison-1.png)

    p_vals <- tibble(
      LinearModel = tidy_mod %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("p.value") %>% unlist,
      ComBat = tidy_mod_combat %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("p.value") %>% unlist,
      SVA = tidy_mod_sva %>% filter(term == "as.factor(bottomly_metadata$strain)DBA/2J") %>% select("p.value") %>% unlist)

    #png("Guzek_problem5b.png", width = 1400, height = 500)
    p_vals.gather <- gather(p_vals)
    ggplot(p_vals.gather, aes(x = value)) + geom_histogram(fill = "darkorange") + facet_wrap(~key)

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](homework_files/figure-markdown_strict/histogram%20os%20pvals%20from%20lm,%20sva,%20and%20combat-1.png)

    #dev.off()
