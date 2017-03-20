#!/usr/bin/env Rscript


elimcor_sansY<-function(X,s=0.95)
    {
	#X matrice contenant les variables à grouper
	#Y vecteur contenant les groupes à prédire
	#s valeur seuil de corrélation
	#print (rownames(X))
  correl=cor(as.matrix(X))
	stop=F
	possetap=1:ncol(X)
	groupes=as.list(1:ncol(X))

	while (stop==F)
	    {
		##regroupement des var pour lesquelles |corr|>0.95
		gplist<-list(NULL)
		possglob=1:ncol(correl)
		for (i in 1:(ncol(correl)))
		   {
			poss=possglob[-i]
			gplist[[i]]=c(i,poss[abs(correl[i,poss])>s])
		   }
		##on trie les groupes du plus gros au plus petit
		gplisteff=unlist(lapply(gplist,length))
		if (any(gplisteff>1))
 		    {
			gplistfin=gplist[gplisteff>1]
			gplistuniq=unlist(gplist[gplisteff==1])
			gpsel=NULL
			##on sélectionne dans chaque groupe une variable au hasard
			for (i in 1:length(gplistfin))
			    {
				selloc=min(gplistfin[[i]])
				gploc=groupes[[possetap[selloc]]]
				for (j in 1:length(gplistfin[[i]]))
				    {
					gploc=c(gploc,groupes[[possetap[gplistfin[[i]][j]]]])				    }
				groupes[[possetap[selloc]]]=unique(gploc)
				gpsel=c(gpsel,selloc)
  			    }
			possetap=possetap[c(gplistuniq,unique(gpsel))]
			correl=cor(X[,possetap])
		    }
		else stop=T	
 	   }
	#groupeseff=unlist(lapply(groupes,length))
	#groupes=groupes[groupeseff>1]
	return(list(possetap=possetap,groupes=groupes))
    }

