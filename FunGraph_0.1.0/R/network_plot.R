#' @title convert ODE results to basic network plot table
#' @param i list result from get_all_net(only use a single list)
#' @return a table with basic information to plot network
get_after <- function(i){
  temp <- matrix(NA,nrow = length(i[[2]]),ncol=3)
  temp[,1] <- i[[2]]
  temp[,2] <- i[[1]]
  temp[,3] <- i[[5]][2:(length(i[[2]])+1)]

  colnames(temp) <- c('from','to','dep_effect')
  temp <- data.frame(temp)
  temp[,3] <- as.numeric(as.character(temp[,3]))
  return(temp)
}

#' @title convert ODE results to basic network plot table for Cytoscape
#' @importFrom utils write.csv
#' @importFrom stats aggregate
#' @param input1 list result from get_ode_par
#' @param input2 list result from get_all_net
#' @param write_data write the result, default is TRUE
#' @return a list with node and link attribute
#' @export
get_net_output <- function(input1, input2, write_data = TRUE){
  data = input1$dataset

  links <- do.call(rbind,lapply(input2, get_after))
  links$effect_type <- sapply(1:nrow(links),function(c)get_effect_type(c,links))

  get_ind <- function(i){
    temp <- i[[5]][1]
    return(temp)
  }

  nodes <- data.frame(unique(links[,2]),row.names(data),sapply(input2,get_ind))
  colnames(nodes) <- c("id","name","received_effect")
  nodes$influence <- aggregate(dep_effect ~ to, data = links, sum)[,2]

  links$dep_effect <- abs(links$dep_effect)
  if (write_data == TRUE) {
    cat('write data for Cytoscape')
    write.csv(links, file = 'links_data.csv')
    write.csv(nodes, file = 'nodes_data.csv')
  }
  return(list(links = links,
              nodes = nodes))
}

#' @title convert >=0 value to '+' and < 0 value = '-'
#' @param i scalar of row number
#' @param links matrix contain links data of network
#' @return a characte 0f '+' or '-'
get_effect_type <- function(i,links){
  tmp <- links$dep_effect[i]
  if (tmp >= 0 ) {
    tmp2 <- '+'
  } else {
    tmp2 <- '-'
  }
  return(tmp2)
}

#' @title calculate maximum effect of a network
#' @importFrom stats aggregate
#' @param k list result from get_all_net
#' @return vector of effect
#' @export
get_max_effect <- function(k){
  after <- do.call(rbind,lapply(k, get_after))
  max_dep_effect <- max(abs(after$dep_effect))

  temp <- aggregate(dep_effect ~ to, data = after, sum)
  all_dep_effect <- max(abs(temp$dep_effect))
  return(c(max_dep_effect,all_dep_effect))
}

#' @title generate color for edge
#' @param i row number
#' @param color_data dataframe of all color for edge
#' @return Hex code
get_colour <- function(i, color_data){
  temp <-  round(i*10)
  temp2 <- adjustcolor(color_data[which(color_data[,2]==temp),1], alpha.f = 0.65)
  return(temp2)
}

#' @title replace cut function results to numeric vector
#' @importFrom stringr str_remove_all str_remove
#' @param i characte of cut function results
#' @return cleaned numeric vector
replace_character <- function(i){
  tmp <- str_remove_all(i, "]")
  tmp1 <- str_remove(tmp, "\\(")
  tmp2 <- as.numeric(unlist(strsplit(tmp1,split='\\,')))
  return(tmp2)
}

#' @title generate color for edge
#' @import igraph
#' @importFrom grDevices adjustcolor colorRampPalette dev.off pdf
#' @param k list result from get_all_net
#' @param title text for plot title
#' @param max_effect matrix for
#' @param color_size scalar control graident size of colour
#' @param save_plot save pdf file
#' @return network plot
#' @export
network_plot <- function(k, title, max_effect, color_size = 100, save_plot = TRUE){
  #ind effect control node size
  get_extra <- function(i){
    temp <- i[[5]][1]
    return(temp)
  }
  extra <- sapply(k,get_extra)
  
  after <- do.call(rbind,lapply(k, get_after))
  
  colfunc <- colorRampPalette(c("#619CFF",
                                "#ffdead",
                                "#F8766D"))
  
  links_col = data.frame(color = colfunc(color_size),
                         y = cut(seq(-max(max_effect[,1]),max(max_effect[,1]),length=color_size),color_size))
  
  links_interval = t(sapply(1:color_size,function(c) replace_character(links_col$y[c])))
  
  links <- after
  colnames(links) <- c("from","to","weight")
  for (i in 1:nrow(links)) {
    links$edge.colour[i] <- links_col$color[which(sapply(1:color_size,function(c)
      findInterval(links$weight[i],c(links_interval[c,])))==1)] #add colour for links
  }
  
  #nodes
  node_col <- data.frame(color = colfunc(color_size),
                         y = cut(seq(-max(max_effect[,2]),max(max_effect[,2]),length=color_size),color_size))
  node_interval = t(sapply(1:color_size,function(c) replace_character(node_col$y[c])))
  nodes <- data.frame(unique(links[,2]),unique(links[,2]),extra)
  colnames(nodes) <- c("id","name","ind_effect")
  nodes$influence <- round(aggregate(weight ~ to, data = links, sum)[,2],3)
  for (i in 1:nrow(nodes)) {
    nodes$colour[i] <- node_col$color[which(sapply(1:color_size,function(c)
      findInterval(nodes$influence[i],c(node_interval[c,])))==1)] #add colour for links
  }
  
  #normalization
  normalization <- function(x){(x-min(x))/(max(x)-min(x))*1.5+0.3}
  
  #final plot
  links[,3] <- normalization(abs(links[,3]))
  nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  net <- graph_from_data_frame( d=links,vertices = nodes,directed = T )
  
  #layout
  l <- layout_randomly(net)
  
  if (save_plot == TRUE) {
    cat('Save PDF plot','\n')
    pdf(paste0(title,"_network_plot.pdf"),width=10,height=10)
    plot.igraph(net,
                vertex.label=V(net)$name,
                vertex.label.color="black",
                vertex.shape="circle",
                vertex.label.cex=V(net)$ind_effect*0.5,
                vertex.size=V(net)$ind_effect*10,
                edge.curved=0.05,
                edge.color=E(net)$edge.colour,
                edge.frame.color=E(net)$edge.colour,
                edge.width=E(net)$weight*4,
                vertex.color=V(net)$colour,
                layout=l,
                main=title,
                margin=c(-.05,-.05,-.05,-.05)
    )
    dev.off()
  }
  else{
    plot.igraph(net,
                vertex.label=V(net)$name,
                vertex.label.color="black",
                vertex.shape="circle",
                vertex.label.cex=V(net)$ind_effect*0.5,
                vertex.size=V(net)$ind_effect*10,
                edge.curved=0.05,
                edge.color=E(net)$edge.colour,
                edge.frame.color=E(net)$edge.colour,
                edge.width=E(net)$weight*4,
                vertex.color=V(net)$colour,
                layout=l,
                main=title,
                margin=c(-.05,-.05,-.05,-.05)
    )
  }
}

