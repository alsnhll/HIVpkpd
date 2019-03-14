single_mut_plot <- function(plot_type,window_units,d,Cwt,twt,tmax,df,axis_units,show_muts){

# Plot Outcomes:
if(window_units == "days"){
  l <- list(
    font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
    orientation = 'v',
    bgcolor = "white",
    bordercolor = "white",
    traceorder = "normal")
  
  tval <- seq(0,tmax,by=0.25)
  ttxt <- rep("",length(tval))  
  ttxt[seq(1,length(ttxt),by=4)] <- as.character(tval)[seq(1,length(ttxt),4)]
  
  
  p = plot_ly(y="Wild Type",x = twt,type = 'bar',orientation = 'h',width = 1100, height = 550,
              marker = list(color = rgb(162, 208, 121,max=255)),
              name = paste('','Suppressive therapy','',sep = '<br>'))
  p=add_trace(p,y="Wild Type",x=0,type='bar',orientation='h',
              marker=list(color=rgb(127, 20, 22,max=255)),
              name=paste('','Only MUT can grow','',sep='<br>'))
  p = add_trace(p,y ="Wild Type",x = 0,type = 'bar',orientation = 'h',
                marker = list(color = 'indianred'),
                name = paste('','Only WT can grow', '',sep='<br>'))
  p = add_trace(p,y = "Wild Type",x = 0,
                marker=list(color=rgb(237, 28, 36,max=255)),
                name=paste('','MUT + WT growth', 'MUT selected',sep='<br>'))
  p = add_trace(p,y = "Wild Type",x = tmax-twt,
                marker = list(color=rgb(202, 216, 224,max=255)),
                name = paste('','MUT + WT growth', 'WT selected',sep='<br>'))
  p=layout(p,xaxis=list(zeroline=F,  ticks = "outside",tickvals = tval,ticktext = ttxt,title='Time since last dose (days)'),yaxis=list(title=''),hovermode="closest",
           title=paste("Selection regimes with", d$drugname,"treatment"),barmode='stack',
           margin=list(l=120,t=82),showlegend=TRUE,legend = l, autosize = F)
  
  for(i in 1:length(plot_type)){
    # Plot Outcomes:
    if(plot_type[i] == 1){
      p=add_trace(p,y=paste(show_muts[i]),x=rev(df$ttop_mut[i]),type='bar',orientation='h',
                  marker = list(color = rgb(162, 208, 121,max=255)),showlegend = FALSE,
                  name='Suppressive therapy')
      p=add_trace(p,y=paste(show_muts[i]),x=twt-rev(df$ttop_mut[i]),type='bar',orientation='h',
                  marker=list(color=rgb(127, 20, 22,max=255)),showlegend = FALSE,
                  name=paste('','Only MUT can grow','',sep='<br>'))
      p=add_trace(p,y=paste(show_muts[i]),x=rev(df$tbot_mut[i])-twt,type='bar',orientation='h',
                  marker=list(color=rgb(237, 28, 36,max=255)),showlegend = FALSE,
                  name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
      p=add_trace(p,y=paste(show_muts[i]),x=tmax-rev(df$tbot_mut[i]),type='bar',orientation='h',
                  marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                  name=paste('','MUT + WT grow', 'WT selected',sep='<br>')) 
      p}else if(plot_type[i] == 2){
        
        p=add_trace(p,y=paste(show_muts[i]),x=twt,type='bar',orientation='h',
                    marker = list(color = rgb(162, 208, 121,max=255)),showlegend = FALSE,
                    name='Suppressive therapy')
        p=add_trace(p,y=paste(show_muts[i]),x=tmax-twt,type='bar',orientation='h',
                    marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                    name=paste('WT Growth Selected',sep='<br>')) 
        p
        
      }else{
        
        p = add_trace(p,y=paste(show_muts[i]),x = twt,type = 'bar',orientation = 'h',
                      marker = list(color = rgb(162, 208, 121,max=255)),showlegend = FALSE,
                      name = 'Suppressive therapy')
        p = add_trace(p,y = paste(show_muts[i]),x = df$tmut[i]-twt,type = 'bar',orientation = 'h',
                      marker = list(color = 'indianred'),showlegend = FALSE,
                      name = paste('','Only WT can grow', '',sep='<br>'))
        p = add_trace(p,y = paste(show_muts[i]),x = df$ttop_int[i]-df$tmut[i],type = 'bar',orientation = 'h',
                      marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                      name = paste('','MUT + WT grow', 'WT selected',sep='<br>'))
        p = add_trace(p,y = paste(show_muts[i]),x = df$tbot_int[i] - df$ttop_int[i],type = 'bar',orientation = 'h',
                      marker=list(color=rgb(237, 28, 36,max=255)),showlegend = FALSE,
                      name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
        p = add_trace(p,y = paste(show_muts[i]),x = tmax-df$tbot_int[i],type = 'bar',orientation = 'h',
                      marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                      name = paste('','MUT + WT grow', 'WT selected',sep='<br>'))
        p
      
      
    }}}else{ # If Selection Window Input is Concentrations:
      Cwt <- if(axis_units == "µM")Cwt else if(axis_units == "µg/ml")Cwt*d$mw/1000 else Cwt*d$mw
      cmax <- if(axis_units == "µM")d$cmax else if(axis_units == "µg/ml")d$cmax*d$mw/1000 else d$cmax*d$mw
      
      # Plot Wild Type (to create plotly)
      l <- list(
        font = list(
          family = "sans-serif",
          size = 12,
          color = "#000"),
        orientation = 'v',
        bgcolor = "white",
        bordercolor = "white",
        traceorder = "reversed")
      
      p=plot_ly(y="Wild Type",x=Cwt,type='bar',orientation='h',width = 1100, height = 550,
                marker = list(color=rgb(202, 216, 224,max=255)),
                name=paste('','MUT + WT growth', 'WT selected',sep='<br>'))
      p = add_trace(p,y = "Wild Type",x = 0,
                    marker=list(color=rgb(237, 28, 36,max=255)),
                    name=paste('','MUT + WT growth', 'MUT selected',sep='<br>'))
      p = add_trace(p,y ="Wild Type",x = 0,type = 'bar',orientation = 'h',
                    marker = list(color = 'indianred'),
                    name = paste('','Only WT can grow', '',sep='<br>'))
      p=add_trace(p,y="Wild Type",x=0,type='bar',orientation='h',
                  marker=list(color=rgb(127, 20, 22,max=255)),
                  name=paste('','Only MUT can grow','',sep='<br>'))
      p=add_trace(p,y="Wild Type",x=cmax-Cwt,type='bar',orientation='h',
                  marker = list(color = rgb(162, 208, 121,max=255)),
                  name=paste('','Suppressive therapy','',sep = '<br>'))
      p=layout(p,xaxis=list(range = c(0,cmax),title=paste0("Drug Concentration (", axis_units,")")),yaxis=list(title=''),hovermode="closest",
               title=paste("Selection regimes with", d$drugname,"treatment"),barmode='stack',
               margin=list(l=120,t=82),showlegend = TRUE,legend = l, autosize = F)
      
      for(i in 1:length(plot_type)){
      if(plot_type[i] == 1){
        C_intersect <- if(axis_units == "µM")df$C_intersect[i] else if(axis_units == "µg/ml")df$C_intersect[i]*d$mw/1000 else df$C_intersect[i]*d$mw
        Cmut <- if(axis_units == "µM")df$Cmut[i] else if(axis_units == "µg/ml")df$Cmut[i]*d$mw/1000 else df$Cmut[i]*d$mw
        
        p=add_trace(p,y=paste(show_muts[i]),x=C_intersect,type='bar',orientation='h',
                    marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                    name=paste('','MUT + WT grow', 'WT selected',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=Cwt-C_intersect,type='bar',orientation='h',
                    marker=list(color=rgb(237, 28, 36,max=255)),showlegend = FALSE,
                    name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=Cmut-Cwt,type='bar',orientation='h',
                    marker=list(color=rgb(127, 20, 22,max=255)),showlegend = FALSE,
                    name=paste('','Only MUT can grow','',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=cmax-Cmut,type='bar',orientation='h',
                    marker = list(color = rgb(162, 208, 121,max=255)),showlegend = FALSE,
                    name='Suppressive therapy')
        p
        
        
      }else{if(plot_type ==2){
        
        p=add_trace(p,y=paste(show_muts[i]),x=Cwt,type='bar',orientation='h',
                    marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                    name=paste('','MUT + WT grow', 'WT selected',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=cmax-Cwt,type='bar',orientation='h',
                    marker = list(color = rgb(162, 208, 121,max=255)),showlegend = FALSE,
                    name='Suppressive therapy')
        p
        
        
      }else{
        Cint_bot <- if(axis_units == "µM")df$Cint_bot[i] else if(axis_units == "µg/ml")df$Cint_bot[i]*d$mw/1000 else df$Cint_bot[i]*d$mw
        Cint_top <- if(axis_units == "µM")df$Cint_top[i] else if(axis_units == "µg/ml")df$Cint_top[i]*d$mw/1000 else df$Cint_top[i]*d$mw
        Cmut <- if(axis_units == "µM")df$Cmut[i] else if(axis_units == "µg/ml")df$Cmut[i]*d$mw/1000 else df$Cmut[i]*d$mw
        
        p=add_trace(p,y=paste(show_muts[i]),x=Cint_bot,type='bar',orientation='h',
                    marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                    name=paste('','MUT + WT grow', 'WT selected',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=Cint_top-Cint_bot,type='bar',orientation='h',
                    marker=list(color=rgb(237, 28, 36,max=255)),showlegend = FALSE,
                    name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=Cmut-Cint_top,type='bar',orientation='h',
                    marker = list(color=rgb(202, 216, 224,max=255)),showlegend = FALSE,
                    name = paste('','MUT + WT grow', 'WT selected',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=Cwt-Cmut,type='bar',orientation='h',
                    marker = list(color = 'indianred'),showlegend = FALSE,
                    name = paste('','Only WT can grow', '',sep='<br>'))
        p=add_trace(p,y=paste(show_muts[i]),x=cmax-Cwt,type='bar',orientation='h',
                    marker = list(color = rgb(162, 208, 121,max=255)),showlegend = FALSE,
                    name='Suppressive therapy')
        p
        
      }
      }
      }
    }
  p
}