## Statistics RT-qPCR

#### Unprocessed alone
{
    datap <- "./Data/ASF1_chaperone/RT-qPCR_data_pre_reshaped.txt"
    dt <- read.table(datap, h=T, stringsAsFactors=F, sep='\t', dec=',')

    dt$target <- factor(dt$target, levels=unique(dt$target))
    dt$siRNA <- factor(dt$siRNA, levels=unique(dt$siRNA))

    dt$log2_foldchange <- log2(dt$FoldChange)


    #### comparison to theoric
    for (target in levels(dt$target)) {
        for (siRNA in levels(dt$siRNA)) {
            print(c(target, siRNA))

            sdt <- dt[dt$target == target & dt$siRNA == siRNA,]
            mstd <- sd(dt$log2_foldchange)

            print(t.test(formula=log2_foldchange ~ 1, data=sdt, alternative="two.sided"))
        }
    }


    #### comparison between conditions
    for (sii in 1:(length(levels(dt$siRNA)) - 1)) {
        for (sij in (sii + 1):length(levels(dt$siRNA))) {
            sia <- levels(dt$siRNA)[sii]
            sib <- levels(dt$siRNA)[sij]

            for (target in levels(dt$target)) {
                print(c(target, sia, sib))

                sdta <- dt[dt$target == target & dt$siRNA == sia,]
                sdtb <- dt[dt$target == target & dt$siRNA == sib,]
                print(t.test(sdta$log2_foldchange, sdtb$log2_foldchange), alternative="two.sided")
            }
        }
    }

}


####
