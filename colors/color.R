get_domain_cols = function(){
  c("#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#f19066",
    "#f5cd79", "#546de5", "#e15f41", "#c44569", "#786fa6", "#f8a5c2",
    "#63cdda", "#ea8685", "#596275", "#574b90", "#f78fb3", "#3dc1d3",
    "#e66767", "#303952")
}

get_vcColors = function(alpha = 1, websafe = FALSE, named = TRUE){
  if(websafe){
    col = c("#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3",
            "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", "#CDDC39",
            "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548", "#9E9E9E",
            "#607D8B")
  }else{
    col = c(RColorBrewer::brewer.pal(11, name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060', '#535c68')
    col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  }

  if(named){
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                          'RNA','Splice_Site','Intron','Frame_Shift_Ins','In_Frame_Del','ITD','In_Frame_Ins',
                          'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event', 'pathway')
  }

  col
}

get_titvCol = function(alpha = 1){
  col = c("#F44336", "#3F51B5", "#2196F3", "#4CAF50", "#FFC107", "#FF9800")
  #col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'forestgreen', 'deeppink3')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  col
}
