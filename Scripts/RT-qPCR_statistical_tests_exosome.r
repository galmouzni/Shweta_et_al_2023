## Statistics RT-qPCR

#### Unprocessed alone
{
    datap <- "./Data/ASF1_chaperone/RT-qPCR_data_exosome_reshaped.txt"
    dt <- read.table(datap, h=T, stringsAsFactors=F, sep='\t', dec=',')

    dt$target <- factor(dt$target, levels=unique(dt$target))
    dt$signal <- factor(dt$signal, levels=unique(dt$signal))
    dt$siASF1 <- factor(dt$siASF1, levels=unique(dt$siASF1))
    dt$si44 <- factor(dt$si44, levels=unique(dt$si44))

    dt$log2_foldchange <- log2(dt$foldchange)


    #### comparison to theoric
    for (targ in levels(dt$target)) {
        for (si44s in levels(dt$si44)) {
            for (sgl in levels(dt$signal)) {
                print(c(targ, si44s, sgl))

                sdt <- dt[dt$target == targ & dt$si44 == si44s & dt$signal == sgl,]
                mstd <- sd(dt$log2_foldchange)

                print(t.test(formula=log2_foldchange ~ 1, data=sdt, alternative="two.sided"))
            }
        }
    }


    #### comparison between conditions
    for (sii in 1:(length(levels(dt$si44)) - 1)) {
        for (sij in (sii + 1):length(levels(dt$si44))) {
            sia <- levels(dt$si44)[sii]
            sib <- levels(dt$si44)[sij]

            for (target in levels(dt$target)) {
                for (sgl in levels(dt$signal)) {
                    print(c(sia, sib, target, sgl))

                    sdta <- dt[dt$target == target & dt$si44 == sia & dt$signal == sgl,]
                    sdtb <- dt[dt$target == target & dt$si44 == sib & dt$signal == sgl,]
                    print(t.test(sdta$log2_foldchange, sdtb$log2_foldchange), alternative="two.sided")
                }
            }
        }
    }

}


####
