#-----R script to import functions I have written used in the MOFA2 analysis
#-----Amir Chaman Baz 06_06_24
#-------------------------------------------------------------------------------
#Imports:   corr_factors_cov_spearman, Corr_factors_R2, plot_weights_mapped,
#           plot_data_scatter_mapped, mapProtein, mapRNA,  makecorplot,  
#           regressor, rankplot, IQR_Filter
#-------------------------------------------------------------------------------
  

  
#----- IQR filtering of clinical outcomes, from: Phenotype_filter.R DvA
#-----This function removes outliers from a vector based on the 
#-----Interquartile Range Rule
#-----The upper and lower boundaries are respectively Q1-n*IQR and Q3+n*IQR
      IQR_filter <- function(x, n){

      Q <- quantile(x, probs=c(0.25, 0.75), na.rm = T)
      iqr <- IQR(x, na.rm = T)
      low <- unname(Q[1] - n*iqr)
      up <- unname(Q[2] + n*iqr)
      x[x < low | x > up] <- NA
      return(x)
      }
  
  
  
#-----Function to correlate factors with covariates (Spearman)
#-----Returns heatmap (ggcorrplot)
corr_factors_cov_spearman<-function (object, 
                                     covariates, 
                                     factors = "all", 
                                     groups = "all", 
                                     abs = FALSE, 
                                     plot = c("r", "pval"), 
                                     alpha = 0.05, 
                                     return_data = FALSE, 
                                     transpose = FALSE, ...)   
{
       
      groups <- MOFA2:::.check_and_get_groups(object, groups)
      metadata <- samples_metadata(object)
      metadata <- metadata[metadata$group %in% groups, ]
      if (is.character(covariates)) {
            stopifnot(all(covariates %in% colnames(metadata)))
            covariates <- metadata[, covariates, drop = FALSE]
      }
      else if (is.data.frame(covariates)) {
            samples <- metadata$sample
            stopifnot(all(rownames(covariates) %in% samples))
            covariates <- metadata[match(rownames(covariates), metadata$sample), 
                                   ]
      }
      cols <- which(sapply(covariates, is.character))
      if (length(cols >= 1)) {
            covariates[cols] <- lapply(covariates[cols], as.factor)
      }
      cols <- which(!sapply(covariates, class) %in% c("numeric", 
                                                      "integer"))
      if (length(cols >= 1)) {
            cols.factor <- which(sapply(covariates, class) == "factor")
            covariates[cols] <- lapply(covariates[cols], as.numeric)
            covariates[cols] <- lapply(covariates[cols], as.numeric)
      }
      stopifnot(all(sapply(covariates, class) %in% c("numeric", 
                                                     "integer")))
      factors <- MOFA2:::.check_and_get_factors(object, factors)
      Z <- get_factors(object, factors = factors, groups = groups, 
                       as.data.frame = FALSE)
      Z <- do.call(rbind, Z)
      cor <- psych::corr.test(Z, covariates, method = "spearman", 
                              adjust = "BH") 
      plot <- match.arg(plot)
      if (plot == "r") {
            stat <- cor$r
            if (abs) 
                  stat <- abs(stat)
            if (transpose) 
                  stat <- t(stat)
            if (return_data) 
                  return(stat)
            return(ggcorrplot(stat, tl.col = "black", title = "Spearman Correlation", 
                       ...))
      }
	  if (plot == "pval") {
            stat <- cor$p.adj
            if (return_data) 
                  return(stat)
            return(ggcorrplot(stat, tl.col = "black", title = "Spearman Correlation",  ...))
      }
}



#-----Function to fit lm to factor
#-----Returns a list(z, covariates) for further processing
Corr_factors_R2<-function (object, 
                           covariates, 
                           factors = "all", 
                           groups = "all", 
                           abs = FALSE, 
                           plot = c("r"), 
                           alpha = 0.05, 
                           return_data = FALSE, 
                           transpose = FALSE, ...) 
{
      groups <- MOFA2:::.check_and_get_groups(object, groups)
      metadata <- samples_metadata(object)
      metadata <- metadata[metadata$group %in% groups, ]
      if (is.character(covariates)) {
            stopifnot(all(covariates %in% colnames(metadata))) 
            covariates <- metadata[, covariates, drop = FALSE]
      }
      else if (is.data.frame(covariates)) {
            samples <- metadata$sample
            stopifnot(all(rownames(covariates) %in% samples))
            covariates <- metadata[match(rownames(covariates), metadata$sample), 
                                   ]
      }
      cols <- which(sapply(covariates, is.character))
      if (length(cols >= 1)) {
            covariates[cols] <- lapply(covariates[cols], as.factor)
      }
      cols <- which(!sapply(covariates, class) %in% c("numeric", 
                                                      "integer"))
      if (length(cols >= 1)) {
            cols.factor <- which(sapply(covariates, class) == "factor")
            covariates[cols] <- lapply(covariates[cols], as.numeric)
            covariates[cols] <- lapply(covariates[cols], as.numeric)
      } 
      stopifnot(all(sapply(covariates, class) %in% c("numeric", 
                                                     "integer")))
      factors <- MOFA2:::.check_and_get_factors(object, factors)
      Z <- get_factors(object, factors = factors, groups = groups, 
                       as.data.frame = FALSE)
      Z <- do.call(rbind, Z)
      return(list(Z, covariates))
      
}

#-----Plot Weights (tables) mapped to either protein or genenames
#-----Returns the plots (ggobject)
plot_weights_mapped<- function (object, 
                                view = 1, 
                                factors = 1, 
                                nfeatures = 10, 
                                abs = TRUE, 
                                scale = TRUE, 
                                sign = "all", 
                                genenames=proteinnames){
      if (!is(object, "MOFA")) 
            stop("'object' has to be an instance of MOFA")
      if (nfeatures <= 0) 
            stop("'nfeatures' has to be greater than 0")
      if (sign == "all") {
            abs <- TRUE
      }
      if (is.numeric(view)) 
            view <- views_names(object)[view]
      stopifnot(view %in% views_names(object))
      view <- MOFA2:::.check_and_get_views(object, view)
      factors <- MOFA2:::.check_and_get_factors(object, factors)
      
      W <- get_weights(object, 
                       factors = factors, 
                       views = view, 
                       as.data.frame = TRUE)
      if (scale) 
            W$value <- W$value/max(abs(W$value))
      W <- W[W$value != 0, ]
      W$sign <- ifelse(W$value > 0, "+", "-")
      if (sign == "positive") {
            W <- W[W$value > 0, ]
      }
      else if (sign == "negative") {
            W <- W[W$value < 0, ]
      }
      if (abs) 
            W$value <- abs(W$value)
      W <- W[with(W, order(-abs(value))), ]
      W <- as.data.frame(top_n(group_by(W, factor), n = nfeatures, 
                               wt = value))
      W$feature_id <- W$feature
      if ((length(unique(W$view)) > 1) && 
          (nfeatures > 0) && 
          (any(duplicated(W[W$factor == factors[1], ]$feature_id)))) {
            message("Duplicated feature names across views, we will add the view name as a prefix")
            W$feature_id <- paste(W$view, W$feature, sep = "_")
      }
      W$feature_id <- factor(W$feature_id, levels = rev(unique(W$feature_id)))
      if(view=="rna"){
            W$feature_id<- as.character(W$feature_id)
            for(i in W$feature_id)
            {
                  W$feature_id[which(W$feature_id==i)]<- genenames$gene_name[which(genenames$gene_id==i)]
            }
            W$feature_id <- factor(W$feature_id, levels = rev(unique(W$feature_id)))
      }
      if(view=="protein"){
            W$feature_id<- as.character(W$feature_id)
            
            for(i in W$feature_id)
            {
                  
                  xx<- strsplit(i, ";")
                  
                  replac<-toString(proteinnames$protein.name[which(proteinnames$protein.id==unlist(xx))])
                  if(!(replac=="")){
                        W$feature_id[which(W$feature_id==i)]<- replac
                  }
            }
            W$feature_id <- factor(W$feature_id, levels = rev(unique(W$feature_id)))
      }
      p <- ggplot(W, aes_string(x = "feature_id", y = "value")) + 
            geom_point(size = 2) + 
            geom_segment(aes_string(xend = "feature_id"), 
                         size = 0.75, yend = 0) + 
            scale_colour_gradient(low = "grey", 
                                  high = "black") + 
            coord_flip() + 
            labs(y = "Weight") + 
            theme_bw() + 
            theme(axis.title.x = element_text(color = "black"), 
                  axis.title.y = element_blank(), 
                  axis.text.y = element_text(size = rel(1.1), 
                                             hjust = 1, 
                                             color = "black"), 
                  axis.text.x = element_text(color = "black"), 
                  axis.ticks.y = element_blank(), axis.ticks.x = element_line(), 
                  legend.position = "top", legend.title = element_blank(), 
                  legend.text = element_text(color = "black"), 
                  legend.key = element_rect(fill = "transparent"), 
                  strip.text = element_text(size = rel(1.2)), 
                  panel.background = element_blank(), 
                  panel.spacing = unit(1, "lines"), 
                  panel.grid.major.y = element_blank(), 
            ) + 
            facet_wrap(~factor, nrow = 1, scales = "free")
      
      if (sign == "negative") 
            p <- p + scale_x_discrete(position = "top")
      if (abs) {
            p <- p + 
                  ylim(0, max(W$value) + 0.1) + 
                  geom_text(label = W$sign, 
                            y = max(W$value) + 0.1, 
                            size = 10)
      }
      return(p)
}


#-----Scatter plots with mapped names (protein, rna) 
#-----Returns ggobject
plot_data_scatter_mapped<- function (object, 
                                     factor = 1, 
                                     view = 1, 
                                     groups = "all", 
                                     features = 10, 
                                     sign = "all", 
                                     color_by = "group", 
                                     legend = TRUE, 
                                     alpha = 1, 
                                     shape_by = NULL, 
                                     stroke = NULL, 
                                     dot_size = 2.5, 
                                     text_size = NULL, 
                                     add_lm = TRUE, 
                                     lm_per_group = TRUE, 
                                     imputed = FALSE, 
                                     genenames=proteinnames) 
{
      if (!is(object, "MOFA")) 
            stop("'object' has to be an instance of MOFA")
      stopifnot(length(factor) == 1)
      stopifnot(length(view) == 1)
      if (lm_per_group) 
            add_lm = TRUE
      groups <- MOFA2:::.check_and_get_groups(object, groups)
      factor <- MOFA2:::.check_and_get_factors(object, factor)
      view <- MOFA2:::.check_and_get_views(object, view)
      N <- get_dimensions(object)[["N"]]
      W <- get_weights(object)[[view]][, factor]
      if (imputed) {
            Y <- do.call(cbind, object@imputed_data[[view]][groups])
      }
      else {
            Y <- do.call(cbind, object@data[[view]][groups])
      }
      Z <- get_factors(object, factors = factor, groups = groups, 
                       as.data.frame = TRUE)
      Z <- Z[, c("sample", "value")]
      colnames(Z) <- c("sample", "x")
      if (sign == "all") {
            W <- abs(W)
      }
      else if (sign == "positive") {
            W <- W[W > 0]
      }
      else if (sign == "negative") {
            W <- W[W < 0]
      }
      if (is(features, "numeric")) {
            if (length(features) == 1) {
                  features <- names(tail(sort(abs(W)), n = features))
            }
            else {
                  features <- names(sort(-abs(W))[features])
            }
            stopifnot(all(features %in% features_names(object)[[view]]))
      }
      else if (is(features, "character")) {
            stopifnot(all(features %in% features_names(object)[[view]]))
      }
      else {
            stop("Features need to be either a numeric or character vector")
      }
      W <- W[features]
      if (length(color_by) == 1 & is.character(color_by)) 
            color_name <- color_by
      if (length(shape_by) == 1 & is.character(shape_by)) 
            shape_name <- shape_by
      color_by <- MOFA2:::.set_colorby(object, color_by)
      shape_by <- MOFA2:::.set_shapeby(object, shape_by)
      df1 <- merge(Z, color_by, by = "sample")
      df1 <- merge(df1, shape_by, by = "sample")
      foo <- list(features)
      names(foo) <- view
      if (isTRUE(imputed)) {
            df2 <- get_imputed_data(object, groups = groups, views = view, 
                                    features = foo, as.data.frame = TRUE)
      }
      else {
            df2 <- get_data(object, groups = groups, features = foo, 
                            as.data.frame = TRUE)
      }
      df2$sample <- as.character(df2$sample)
      df <- dplyr::left_join(df1, df2, by = "sample")
      df <- df[!is.na(df$value), ]
      if(view=="rna"){
            df$feature<- as.character(df$feature)
            for(i in df$feature)
            {
                  df$feature[which(df$feature==i)]<- genenames$gene_name[which(genenames$gene_id==i)]
            }
            df$feature <- factor(df$feature, levels = rev(unique(df$feature)))
      }
      if(view=="protein"){
            
            df$feature<- as.character(df$feature)
            
            for(i in df$feature)
            {
                  xx<- strsplit(i, ";")
                  replac<-toString(proteinnames$protein.name[which(proteinnames$protein.id==unlist(xx))])
                  if(!(replac=="")){
                        df$feature[which(df$feature==i)]<- replac
                  }
            }
            df$feature <- factor(df$feature, levels = rev(unique(df$feature)))
      }
      
      if (is.null(stroke)) {
            stroke <- MOFA2:::.select_stroke(N = length(unique(df$sample)))
      }
      if (add_lm && is.null(text_size)) {
            text_size <- MOFA2:::.select_pearson_text_size(N = length(unique(df$feature)))
      }
      axis.text.size <- MOFA2:::.select_axis.text.size(N = length(unique(df$feature)))
      p <- ggplot(df, aes_string(x = "x", y = "value")) + 
            geom_point(aes_string(fill = "color_by", shape = "shape_by"), 
                       colour = "black", 
                       size = dot_size, 
                       stroke = stroke, 
                       alpha = alpha) + 
            labs(x = "Factor values", y = "") + 
            facet_wrap(~feature, scales = "free_y") + 
            theme_classic() + 
            theme(axis.text = element_text(size = rel(axis.text.size), 
                                                             color = "black"), 
                  axis.title = element_text(size = rel(1), color = "black"), 
                  strip.background = element_blank())      
            if (add_lm) {
                  if (lm_per_group && length(groups) > 1) {
                  p <- p + stat_smooth(formula = y ~ x, aes_string(color = "group"), 
                                       method = "lm", alpha = 0.4) + 
                        ggpubr::stat_cor(aes_string(color = "group", 
                                                    label = "..r.label.."), 
                                         method = "pearson",
                                         label.sep = "\n", 
                                         output.type = "latex", 
                                         size = text_size)
            }
            else {
                  p <- p + stat_smooth(formula = y ~ x, method = "lm", 
                                       color = "grey", fill = "grey", alpha = 0.4) #+ 
                  #ggpubr::stat_cor(method = "pearson", label.sep = "\n", 
                  #output.type = "latex", size = text_size, color = "red")
            }
      }
      p <- MOFA2:::.add_legend(p, df, legend, color_name, shape_name)
      return(p)
}


#-----Function to map protein Id's to protein names. 
#-----When a protein identificator consists of several id's, map the single id's 
#-----and return a concatenated string of names separated by ";". 
#-----If no name is found, the id original is returned again
mapProtein<- function(proteinlist, proteinname=proteinnames){
      list1<- list()     
      for(i in proteinlist){
            if(length(grep(";", i))>0){
                  protname<-paste(unlist(lapply(unlist(strsplit(i, ";")), 
                        function(x){
                              matchname<-proteinnames$protein.name[which(proteinnames$protein.id==x)]
                        if(length(matchname)>0){
                              matchname
                        }else{
                              x
                        }
                                                      }
                        )
                        ), collapse=";")
                  if(length(protname)>0)
                  {
                        list1<-c(list1, protname)
                  }else{
                        list1<-c(list1, i)
                  }
                  
            }else{
                  protname<-proteinname$protein.name[which(proteinname$protein.id==i)]
                  
                  if(length(protname>0)){
                        list1<- c(list1, protname)
                  }else{
                        list1<-c(list1, i)
                  }
                  
            }
            
      }
      return(unlist(list1))
}

mapRna<- function(rnalist, genename=genenames){
      list1<- list()     
      for(i in rnalist){
                  rnaname<-genename$gene_name[which(genename$gene_id==i)]
                  
                  if(length(rnaname>0)){
                        list1<- c(list1, rnaname)
                  }else{
                        list1<-c(list1, i)
                  }
                  
            }
            
      
      return(unlist(list1))
}

#----------Makes correlation plots between feature ~outcome  
#----------Returns single ggplot
makecorplot<-function(sfeature, soutcome, colby="CTGDiagnostic", returncorr=FALSE){
    
    selectedfeature<-rpc2[which(rownames(rpc2)==sfeature),]
    selectedoutcome<-metadatadf[which(rownames(metadatadf)==soutcome),]
    colfill<-metadatadf[which(rownames(metadatadf)==colby),]
    
    selectedoutcome<-selectedoutcome[,colnames(selectedoutcome) %in% 
                                         colnames(selectedfeature)]
    colfill<- colfill[, colnames(colfill) %in% colnames(selectedoutcome)]
    regress<- rbind(selectedoutcome, selectedfeature, colfill)
    regress<- t(regress)
    regress<- as.data.frame(regress)
    
    
    if(!returncorr){
        colnames(regress)<- c("x", "y", "colfill")
    }else{
        colnames(regress)<- c("x", sfeature, "colfill")
    }
    regress$x<- as.numeric(regress$x)
    if(!returncorr){
        regress$y<- as.numeric(regress$y)
    }
    regress$colfill<- as.numeric(regress$colfill)
    
    #regressmodel<- lm(y ~  x, data=regress)
    if(returncorr){ 
        p<-regress
    }else{
        if(soutcome %in% rownames(clinoutcomes_scaled)){
            
            p<- ggplot(regress, aes(x=x, y=y)) +
                geom_point(aes(col=as.factor(colfill)))+
                geom_smooth(method="lm", se=FALSE, color="347FC4")+
                theme_minimal() +
                ggpubr::stat_cor(aes_string( 
                    label = "..r.label.."), method = "spearman", 
                    label.sep = "\n", output.type = "latex", size = 4, color="#61D04F", position = "jitter", hjust=-2.8, vjust=2 )+
                theme(#panel.border=element_rect(colour="black", fill=NA, size=2),
				text=element_text(family = "baskerville") )+
                labs(x="", y="") +
                scale_x_continuous(limits=symmetric_limits)+
                scale_y_continuous(limits=symmetric_limits)+
                scale_color_manual(values=c("#440154FF","#347FC4"))+
                ggtitle(paste0(mapProtein(sfeature))) +
                labs(col=colby) 
        }else{ 
            p<- ggplot(regress, aes(x=x, y=y)) +
                geom_point(aes(col=as.factor(colfill)))+
                geom_smooth(method="lm", se=FALSE, color="347FC4")+
                ggpubr::stat_cor(aes_string( 
                    label = "..r.label.."), method = "spearman", 
                    label.sep = "\n", output.type = "latex", size = 4, color="#61D04F", position = "jitter", hjust=-2.8, vjust=2)+
                theme_minimal() +
                theme(#panel.border=element_rect(colour="black", fill=NA, size=2),
				text=element_text(family = "baskerville"))+
                labs(x="", y="") +
                scale_color_manual(values=c("#440154FF","#347FC4"))+
                ggtitle(paste0(mapProtein(sfeature)))+
                labs(col=colby)
            
        }
    }
    return(p)
}
 
#--------Retrieves features with top weights and maps these to 
#--------function makecorplot, arranges ggobjects in one figure
regressor<- function(soutcome, nfactor, nfeature=4, colby, returncorr=FALSE, nrown=2, ncoln=2){
      
            topweights_protein<- get_weights(MOFAobject_ran, 
                                             views="protein", 
                                             factors=nfactor)
            colnames(topweights_protein$protein)<- "factor"
            topweights_protein<- as.data.frame(topweights_protein)
            topweights_protein$protein<- rownames(topweights_protein)
            topweights_protein$factor<- abs(topweights_protein$factor)
            topweights_protein<- topweights_protein[order(topweights_protein$factor, decreasing = TRUE),]
			
			if(typeof(nfeature)=='double'){
			t<-map2(topweights_protein$protein[1:nfeature],  
            soutcome, colby, returncorr, .f = makecorplot) 
			}else{
					namegene<- nfeature
					nfeature<- length(namegene)
					print(namegene)
					t<-map2(namegene, 
									soutcome, 
									colby, returncorr,
									.f = makecorplot)
    } 
	
            
            if(returncorr){
			corrframe<- data.frame(x=t[[1]]$x)
			for(i in 1:length(t))
			{
				ycol<- as.data.frame(as.numeric(t[[i]][,2]))			
				colnames(ycol)<-colnames(t[[i]])[2]
				corrframe<- cbind(corrframe, ycol)
			}
			
			xdf<- data.frame(x=corrframe$x)
			colnames(xdf)<- soutcome
			colnames(corrframe)<- mapProtein(colnames(corrframe))
			pp<-corr.test(x=xdf, y=corrframe[,-1], method="spearman", adjust="BH")
			return(pp)
			}else{
            #nrown<- floor(sqrt(nfeature))
			#ncoln<- ceiling(nfeature/nrown)
			#if(!(nrown*nrown>=nfeature)){nrown<- nrown+1}
            pp<-ggarrange(plotlist = t, 
                          nrow=nrown, 
                          ncol=ncoln, 
                          common.legend = TRUE)
            if(soutcome=="SMWT"){soutcome<-"6MWT"}
			if(soutcome=="V2Mode"){soutcome<-"CTG-repeat length"}
            pp<-grid.arrange(pp, 
                             left=textGrob(paste0("Protein Expression (Log2)"),
                                           rot=90, gp=gpar(fontfamily="baskerville", 
                                                           fontsize=12)), 
                             bottom=textGrob(paste0(soutcome, "(Z-score)"), 
                                             gp=gpar(fontfamily="baskerville", 
                                                     fontsize=12)))
      
            return(pp)
			} 
} 
  
#-----Function to map rank lists(feature leads)
#-----Returns plot
rankplot<- function(a,b, labels.offset=0.1, arrow.len=0.1)
      {
      old.par<- par(mar=c(1,1,1,1))
      len.1<- length(a)
      len.2<- length(b) 
      
      plot(rep(1, len.1), 1:len.1, pch=20, cex=0.8,
           xlim=c(0,3), ylim=c(max(len.1, len.2),0),
           axes=F, xlab="", ylab="")
      points(rep(1,len.2), 1:len.2, pch=20, cex=0.8)
      
      text(rep(0.94, len.1), len.1:1, len.1:1)
      text(rep(2+labels.offset, len.2), 1:len.2, b, pos = 4)
      
      a.to.b<- match(a,b)
      arrows(rep(1.02, len.1), 
             1:len.1, 
             rep(1.98, len.2), 
             a.to.b, 
             length=arrow.len, 
             angle=40, 
             col = "blue")
      par(old.par)
      
      }

	  
	  
	  
	  
makecorplot_rna<-function(sfeature, soutcome, colby="timepoint", addr=FALSE, returncorr=FALSE){
      
      selectedfeature<-rna_counts2[which(rownames(rna_counts2)==sfeature),]
      selectedoutcome<-metadatadf[which(rownames(metadatadf)==soutcome),]
      colfill<-metadatadf[which(rownames(metadatadf)==colby),]
      
      selectedoutcome<-selectedoutcome[,colnames(selectedoutcome) %in% 
                                             colnames(selectedfeature)]
      colfill<- colfill[, colnames(colfill) %in% colnames(selectedoutcome)]
      regress<- rbind(selectedoutcome, selectedfeature, colfill)
      regress<- t(regress)
      regress<- as.data.frame(regress)
      
      
      if(!returncorr){
      colnames(regress)<- c("x", "y", "colfill")
	  }else{
	  colnames(regress)<- c("x", sfeature, "colfill")
	  }
      regress$x<- as.numeric(regress$x)
	  if(!returncorr){
      regress$y<- as.numeric(regress$y)
	  }
      regress$colfill<- as.numeric(regress$colfill)
	  
      #regressmodel<- lm(y ~  x, data=regress)
       if(returncorr){
	  p<-regress
	  }else{
      if(soutcome %in% rownames(clinoutcomes_scaled)){
            p<- ggplot(regress, aes(x=x, y=y)) +
                  geom_point(aes(col=as.factor(colfill)))+
                  geom_smooth(method="lm", se=FALSE, color="347FC4")+
                  ggpubr::stat_cor(aes_string( 
                                label = "..r.label.."), method = "spearman", 
                     label.sep = "\n", output.type = "latex", size = 4, color="#61D04F", position = "jitter", hjust=-2.8, vjust=2)+
				  theme_minimal() +
                  theme(text=element_text(family = "baskerville") 
                        )+
                  labs(x="", y="") +
                scale_color_manual(values=c("#440154FF","#347FC4"))+
                  scale_x_continuous(limits=symmetric_limits)+
                  scale_y_continuous(limits=symmetric_limits)+
                  ggtitle(paste0(mapRna(sfeature))) +
                  labs(col=colby)
      }else{
            p<- ggplot(regress, aes(x=x, y=y)) +
                  geom_point(aes(col=as.factor(colfill)))+
                  geom_smooth(method="lm", se=FALSE, color="347FC4")+
				  ggpubr::stat_cor(aes_string( 
                                label = "..r.label.."), method = "spearman", 
                     label.sep = "\n", output.type = "latex", size = 4, color="#61D04F", position = "jitter", hjust=-2.8, vjust=2)+
                  theme_minimal() +
                  theme(text=element_text(family = "baskerville"))+
                  labs(x="", y="") +
                scale_color_manual(values=c("#440154FF","#347FC4"))+
                  ggtitle(paste0(mapRna(sfeature)))+
                  labs(col=colby)
            
      }
			}
      return(p)
}
 
regressor_rna<- function(soutcome, nfactor, nfeature=4, colby, returncorr=FALSE, ncoln=2, nrown=2){
    
    topweights_rna<- get_weights(MOFAobject_ran, 
                                     views="rna", 
                                     factors=nfactor)
    colnames(topweights_rna$rna)<- "factor"
    topweights_rna<- as.data.frame(topweights_rna)
    topweights_rna$rna<- rownames(topweights_rna)
    topweights_rna$factor<- abs(topweights_rna$factor)
    topweights_rna<- topweights_rna[order(topweights_rna$factor, decreasing = TRUE),]
    
    if(typeof(nfeature)=='double'){
    t<-map2(topweights_rna$rna[1:nfeature], 
            soutcome, colby, addr, returncorr, .f = makecorplot_rna) 
    }else{
      namegene<- nfeature
      nfeature<- length(namegene)
      t<-map2(namegene, 
            soutcome, colby, addr, returncorr, .f = makecorplot_rna)
    }
	
	if(returncorr){
	
			corrframe<- data.frame(x=t[[1]]$x)
			for(i in 1:length(t))
			{
				ycol<- as.data.frame(as.numeric(t[[i]][,2]))			
				colnames(ycol)<-colnames(t[[i]])[2]
				corrframe<- cbind(corrframe, ycol)
			}
			
			xdf<- data.frame(x=corrframe$x)
			colnames(xdf)<- soutcome
			colnames(corrframe)<- mapRna(colnames(corrframe))
			pp<-corr.test(x=xdf, y=corrframe[,-1], method="spearman", adjust="BH")
			return(pp)
			}else{
    
	#nrown<- floor(sqrt(nfeature))
	#ncoln<- ceiling(nfeature/nrown)
	#if(!(nrown*nrown>=nfeature)){nrown<- nrown+1}
    pp<-ggarrange(plotlist = t, 
                  nrow=nrown, 
                  ncol=ncoln, 
                  common.legend = TRUE)
    if(soutcome=="SMWT"){soutcome<-"6MWT"}
	if(soutcome=="V2Mode"){soutcome<-"CTG-repeat Length"}
    pp<-grid.arrange(pp, 
                     left=textGrob(paste0("RNA Expression (log2)"),
                                   rot=90, gp=gpar(fontfamily="baskerville", 
                                                   fontsize=12)), 
                     bottom=textGrob(paste0(soutcome, "(Z-score)"), 
                                     gp=gpar(fontfamily="baskerville", 
                                             fontsize=12)))
    
    return(pp)
    }
} 


#Further analysis of factor (factorp) hits by pca or pls. functions projects lead hits (nfeatures) and clinical outcomes (noutcomes) 
#in latent space, using either principal component analysis (pca) or partial least squares (pls)
#Maps genes identified in pathways (gProfiler) when view "rna" is used.
#Experimental, not used in report/paper. 
pca_factor<- function(factorp, view="rna", type="pca", nfeatures=20, noutcomes=10, showlabels=TRUE){
  
f11<- get_weights(MOFAobject_ran, views = "all", factors = factorp)
f12<-f11

if(view=="rna"){
f12<-as.data.frame(abs(f12$rna))
f12$feature<- rownames(f12)
colnames(f12)<- c("factor", "feature")
f12<- f12[order(f12$factor, decreasing = TRUE),]

selectf<-rna_counts2[0,]
for(i in f12$feature[1:nfeatures]){
    selectf<- rbind(selectf,rna_counts2[which(rownames(rna_counts2)==i),])
}

selectf<-selectf[, colSums(is.na(selectf))< nrow(selectf)]
#rownames(selectf)<- mapRna(rownames(selectf))

}
if(view=="protein"){
f12<-as.data.frame(abs(f12$protein))
f12$feature<- rownames(f12)
colnames(f12)<- c("factor", "feature")
f12<- f12[order(f12$factor, decreasing = TRUE),]

selectf<-rpc2[0,]
for(i in f12$feature[1:nfeatures]){
    selectf<- rbind(selectf,rpc2[which(rownames(rpc2)==i),])
}

selectf<-selectf[, colSums(is.na(selectf))< nrow(selectf)]
#rownames(selectf)<- mapProtein(rownames(selectf))

}
 
f11_phenotype<-as.data.frame(abs(f11$phenotype))
f11_phenotype$feature<- rownames(f11_phenotype)
colnames(f11_phenotype)<-c("factor", "feature")
f11_phenotype<- f11_phenotype[order(f11_phenotype$factor, decreasing = TRUE),]

phenotypeselect<- clinoutcomes_scaled[0,]
for(i in f11_phenotype$feature[1:noutcomes]){
    phenotypeselect<- rbind(phenotypeselect,clinoutcomes_scaled[which(rownames(clinoutcomes_scaled)==i),])
}
phenotypeselect<- phenotypeselect[, colnames(phenotypeselect) %in% colnames(selectf)]

rnaphenoselect<-rbind(selectf, phenotypeselect)

rnaphenoselect<- t(rnaphenoselect)

if(type=="pls"){
pls.rnapheno<-pls(X=t(selectf), Y=t(phenotypeselect), ncomp = 5)
plsdata<-plotVar(pls.rnapheno, plot=FALSE)
plsdata$pathway<- 0

if(length(gsearesults[[factorp]]$result$intersection)>0 && view=="rna"){
for( i in 1:length(gsearesults[[factorp]]$result$intersection)){
plsdata$pathway[rownames(plsdata) %in% unlist(strsplit(gsearesults[[factorp]]$result$intersection[i], ","))]  <- i
}
}

plsdata$pathway[which(plsdata$Block=="Y")]  <- "clinicaloutcome"

#paste clinicaloutcome names in pathway column
for(i in (nrow(plsdata)-(noutcomes-1)):nrow(plsdata))
{
    plsdata$pathway[i]  <- rownames(plsdata[i,])
}

if(view=="rna"){plsdata$names<- mapRna(plsdata$names)}
if(view=="protein"){plsdata$names<- mapProtein(plsdata$names)}

if(showlabels){
p<- ggplot(plsdata, aes(x=x, y=y, label=names))+
    geom_point(aes(col=pathway))+
    geom_text(aes(col=pathway, hjust=0, vjust=0))
}else{
p<- ggplot(plsdata, aes(x=x, y=y))+
    geom_point(aes(col=pathway))  
}
return(p)
}

if(type=="pca"){
pca.rnapheno<-pca(rnaphenoselect, ncomp = 5)
pcadata<- plotVar(pca.rnapheno, plot=FALSE)
pcadata$pathway<- 0

 
if(length(gsearesults[[factorp]]$result$intersection)>0 && view=="rna"){
for( i in 1:length(gsearesults[[factorp]]$result$intersection)){
pcadata$pathway[rownames(pcadata) %in% unlist(strsplit(gsearesults[[factorp]]$result$intersection[i], ","))]  <- i
}
}

#paste clinicaloutcome names in pathway column
for(i in (nrow(pcadata)-(noutcomes-1)):nrow(pcadata))
{
    pcadata$pathway[i]  <- rownames(pcadata[i,])
}

if(view=="rna"){pcadata$names<- mapRna(pcadata$names)}
if(view=="protein"){pcadata$names<- mapProtein(pcadata$names)}

if(showlabels){
p<- ggplot(pcadata, aes(x=x, y=y, label=names))+
    geom_point(aes(col=pathway))+
    geom_text(aes(col=pathway), hjust=0, vjust=0)
}else{
p<- ggplot(pcadata, aes(x=x, y=y))+
    geom_point(aes(col=pathway))  
}
return(p)
}

}
 
#Reverse matrix for heatmap in figure 1 (no row clustering)
plot_data_heatmap2<- function (object, factor, view = 1, groups = "all", features = 50, 
          annotation_features = NULL, annotation_samples = NULL, transpose = FALSE, 
          imputed = FALSE, denoise = FALSE, max.value = NULL, min.value = NULL, 
          ...) 
{
    if (!is(object, "MOFA")) 
        stop("'object' has to be an instance of MOFA")
    stopifnot(length(factor) == 1)
    stopifnot(length(view) == 1)
    groups <- MOFA2:::.check_and_get_groups(object, groups)
    factor <- MOFA2:::.check_and_get_factors(object, factor)
    view <- MOFA2:::.check_and_get_views(object, view)
    W <- do.call(rbind, get_weights(object, views = view, factors = factor, 
                                    as.data.frame = FALSE))
    Z <- lapply(get_factors(object)[groups], function(z) as.matrix(z[, 
                                                                     factor]))
    Z <- do.call(rbind, Z)[, 1]
    Z <- Z[!is.na(Z)]
    if (isTRUE(denoise)) {
        data <- predict(object, views = view, groups = groups)[[1]]
    }
    else {
        if (isTRUE(imputed)) {
            data <- get_imputed_data(object, view, groups)[[1]]
        }
        else {
            data <- get_data(object, views = view, groups = groups)[[1]]
        }
    }
    if (is(data, "list")) {
        data <- do.call(cbind, data)
    }
    if (is(features, "numeric")) {
        if (length(features) == 1) {
            features <- rownames(W)[tail(order(abs(W)), n = features)]
        }
        else {
            features <- rownames(W)[order(abs(W))[features]]
        }
        features <- names(W[features, ])[order(W[features, ])]
    }
    else if (is(features, "character")) {
        stopifnot(all(features %in% features_names(object)[[view]]))
    }
    else {
        stop("Features need to be either a numeric or character vector")
    }
    data <- data[features, ]
    data <- data[, names(Z)]
    data <- data[, apply(data, 2, function(x) !all(is.na(x)))]
    order_samples <- names(sort(Z, decreasing = TRUE))
    order_samples <- order_samples[order_samples %in% colnames(data)]
    data <- data[, order_samples]
    if (!is.null(annotation_samples)) {
        if (is.data.frame(annotation_samples)) {
            message("'annotation_samples' provided as a data.frame, please make sure that the rownames match the sample names")
            if (any(!colnames(data) %in% rownames(annotation_samples))) {
                stop("There are rownames in annotation_samples that do not correspond to sample names in the model")
            }
            annotation_samples <- annotation_samples[colnames(data), 
                                                     , drop = FALSE]
        }
        else if (is.character(annotation_samples)) {
            stopifnot(annotation_samples %in% colnames(object@samples_metadata))
            tmp <- object@samples_metadata
            rownames(tmp) <- tmp$sample
            tmp$sample <- NULL
            tmp <- tmp[order_samples, , drop = FALSE]
            annotation_samples <- tmp[, annotation_samples, 
                                      drop = FALSE]
            rownames(annotation_samples) <- rownames(tmp)
        }
        else {
            stop("Input format for 'annotation_samples' not recognised ")
        }
        foo <- sapply(annotation_samples, function(x) is.logical(x) || 
                          is.character(x))
        if (any(foo)) 
            annotation_samples[, which(foo)] <- lapply(annotation_samples[, 
                                                                          which(foo), drop = FALSE], as.factor)
    }
    if (!is.null(annotation_features)) {
        stop("'annotation_features' is currently not implemented")
    }
    if (transpose) {
        data <- t(data)
        if (!is.null(annotation_samples)) {
            annotation_features <- annotation_samples
            annotation_samples <- NULL
        }
        if (!is.null(annotation_features)) {
            annotation_samples <- annotation_features
            annotation_features <- NULL
        }
    }
    if (!is.null(max.value)) 
        data[data >= max.value] <- max.value
    if (!is.null(min.value)) 
        data[data <= min.value] <- min.value
    pheatmap::pheatmap(data[nrow(data):1,], annotation_row = annotation_features, annotation_col = annotation_samples, 
             ...)
}