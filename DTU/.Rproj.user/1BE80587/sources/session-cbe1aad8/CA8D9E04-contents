plotDTU=function(df, gene)
  {ggplot(df[df$HGNC==gene,], aes(y = value, x = name, fill = as.character(featureID))) +
    geom_flow(aes(alluvium = featureID), alpha= .5, color = "white",
              curve_type = "linear", 
              width = .5) +
    geom_col(width = .5, color = "white") +
    ggtitle(gene)+scale_fill_brewer(palette="Blues", type="seq", aesthetics=c("colour", "fill")) +ylab("Transcript Usage (%)")+
    theme_bw()}


plotLFC=function(df, gene)
{ggplot(df[df$hgnc==gene,], 
        aes( x=level, y=log2FoldChange, 
             fill=transcript, colour=transcript))+
    geom_col(position="dodge", width=0.5)+
    ggtitle(gene)+
    scale_x_discrete(labels=c("Gene", "Transcript"))+ 
    scale_y_continuous(limits=c(-4,4 ))+
    scale_fill_brewer(palette="Reds", type="seq", aesthetics=c("colour", "fill"))+ 
    ylab("logFC") +
    theme_bw()}

plotDTUcomp=function(df, gene)
  {labels=c("cDNA", "direct RNA")
  ggplot(df[df$HGNC==gene,], aes(y = value, x = name, fill = as.character(featureID))) +
    geom_flow(aes(alluvium = featureID), alpha= .5, color = "grey50",
              curve_type = "linear", 
              width = .5) +
    geom_col(width = .5, color = "grey50") +
    ggtitle(gene)+
    facet_wrap(~protocol)+
    scale_fill_viridis(discrete=TRUE) +
    ylab("Transcript Usage (%)")+ xlab("Procotol")+
    theme_bw()+theme(
      legend.text = element_text(size =14),
      axis.text = element_text(size = 18), 
      axis.title.x=element_text(size=18),
      axis.title.y=element_text(size=18))}
    
plotlengthLFC=function(df, title)
{ggplot(df,aes(y=lfc, x=transcript_length))+ geom_point()+
    ylab("log2 Fold Change")+
      xlab("Transcript Length (bp)")+
    ggtitle(title)+
        facet_wrap(~library)}

spread=function(df)
{newdf=df %>% group_by(groupID) %>% mutate(spread=max(transcript_length)-min(transcript_length))}

plotcounts=function(df, title)
{ggplot(df, aes(x=`cDNA Count`, y=`dRNA Count`))+
    geom_point()+facet_wrap(~quartile)+
    scale_x_continuous(limits=c(0,5000))+
    scale_y_continuous(limits=c(0,5000))+
    geom_smooth(method="lm", se=FALSE)+
    stat_cor(colour="red", geom="text")}

 plotvolcano.DGE=function(df, title)
   {ggplot(df)+ geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig)) + 
     scale_colour_manual(values=volcano.colours) + 
     geom_vline(xintercept=c(-1.5,1.5), colour="grey30") + 
     geom_hline(yintercept=-log10(0.05), colour="grey30") +
     scale_y_continuous(limits=c(0, 10))+
     ggtitle(title)}
   
fgsubset=function(df)
{vec=df[,4]
names(vec)=rownames(df)
return(vec)}

rowsums=function(df) {
  require(dplyr)
  y <- select_if(is_numeric, df)
  rowSums(y, na.rm=T)
}

plot_fc_bambu=function(df, title)
{ggplot(df)+geom_point(aes(x=LFC_FC, y=LFC_bambu, colour=sig))+
    scale_colour_manual(values=sig.colours)+
    geom_vline(xintercept=c(-1.5,1.5), colour="grey30")+
    geom_hline(yintercept=c(-1.5,1.5), colour="grey30")+
    ggtitle(title)
  
  }

