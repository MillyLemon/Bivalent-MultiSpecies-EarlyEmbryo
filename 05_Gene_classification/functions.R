#############################################################################
#Description: Get the coverage from a list of or just one GRanges object
#Notes: Most of the time this is just coverage(src) 
#Author: Matthew Young
#Date Modified: 23/8/2010

getCoverage = function(src){
	if(class(src)=="GRanges"){
		cover=coverage(src)
	}else{
		cover=list()
		for(i in 1:length(src)){
			cover[[i]]=getCoverage(src[[i]])
			message("Got coverage for ",names(src)[i])
		}
		names(cover)=names(src)
	}
	cover
}

#############################################################################
#Description: Extracts the coverage track from a region.  If cover is a list, extracts this region
#Notes: If the chromosome doesn't much, a zero coverage track is returned (think about changing this)
#Author: Matthew Young
#Date Modified: 23/8/2010

getRegion=function(cover,chr,start,end){
	#This is a bit crummy, but will do for now
	if(class(cover)=="SimpleRleList"){
		if(!chr%in%names(cover) | length(cover[[chr]])==0){
			out=Rle(rep(0,end-start+1))
		}else{
			start=min(length(cover[[chr]]),max(1,start))
			end=min(length(cover[[chr]]),max(1,end))
	#		out=seqselect(cover[[chr]],start=start,end=end)
	#		out=window(cover[[chr]],start=start,end=end)
			out=cover[[chr]][IRanges(start=start,end=end)]
		}
	}else{
		out=list()
		for(i in 1:length(cover)){
			out[[i]]=getRegion(cover[[i]],chr=chr,start=start,end=end)
		}
		names(out)=names(cover)
	}
	out
}


#############################################################################
#Description: Bin coverage across the genome using this function
#Notes:
#Author: Matthew Young
#Date Modified: 23/8/2010

binGenome=function(cover,binsize=1000){
	if(class(cover)=="SimpleRleList"){
		tmp=lapply(cover,runsum,k=binsize,endrule="drop")
		out=RleList(lapply(tmp,function(u){u[seq(1,length(u),binsize)]}))
	}else{
		out=list()
		for(i in 1:length(cover)){
			out[[i]]=binGenome(cover[[i]],binsize=binsize)
		}
		names(out)=names(cover)
	}
	out
}


#############################################################################
#Description: Given a coverage track (with the coverage being a pileup over dots), generates a GRanges object with "reads"
#Notes:
#Author: Matthew Young
#Date Modified: 23/8/2010

coverageToReads=function(coverage){
	chrs=rep(names(coverage),sapply(coverage,sum))
	starts=unlist(lapply(coverage,function(u){rep(which(u!=0),as.numeric(u[which(u!=0)]))}))
	GRanges(seqnames=chrs,range=IRanges(start=starts,width=1),strand='*')
}
#############################################################################
#Description: Clusters a set of genes by the data provided
#Notes: Returns the finessed data and the associated metadata
#Data is either a coverage object or a scaled gene object
#Author: Matthew Young
#Date Modified: 28/10/2010

clusterGenes=function(data,up,down,annot,ncluster=5,smoothing=1,cap_quant=.97,cluster=TRUE,plot=TRUE,plot_key=TRUE,status_vector=NULL,exp_vector=NULL,data_vector=NULL,sort_by=c("cluster"),title="Coverage cluster",xlab="% gene length from TSS",xlab_tss="Bases from TSS",status_cols=0:length(unique(status_vector))){
	output=list()
	xlab=xlab
	if(ncluster<2)
		cluster=FALSE
	#the output of cols=colorpanel(64,"blue","white","red").  colorpanel is from the gplots package
	cols=c("#0000FF","#0808FF","#1010FF","#1919FF","#2121FF","#2929FF","#3131FF","#3A3AFF","#4242FF","#4A4AFF","#5252FF","#5A5AFF","#6363FF","#6B6BFF","#7373FF","#7B7BFF","#8484FF","#8C8CFF","#9494FF","#9C9CFF","#A5A5FF","#ADADFF","#B5B5FF","#BDBDFF","#C5C5FF","#CECEFF","#D6D6FF","#DEDEFF","#E6E6FF","#EFEFFF","#F7F7FF","#FFFFFF","#FFFFFF","#FFF7F7","#FFEFEF","#FFE6E6","#FFDEDE","#FFD6D6","#FFCECE","#FFC5C5","#FFBDBD","#FFB5B5","#FFADAD","#FFA5A5","#FF9C9C","#FF9494","#FF8C8C","#FF8484","#FF7B7B","#FF7373","#FF6B6B","#FF6363","#FF5A5A","#FF5252","#FF4A4A","#FF4242","#FF3A3A","#FF3131","#FF2929","#FF2121","#FF1919","#FF1010","#FF0808","#FF0000")
	#Do we have a scaled gene object or a coverage object, if it's a coverage object, we need to know where the TSS is
	message("Fetching/formatting data...")
	if(!missing(annot)){
		xlab=xlab_tss
		if(ncol(values(annot))>0){
			row.ids=values(annot)[,1]
		}else{
			row.ids=1:length(annot)
		}
		#Data fetching is only supported for TSS centric
		range=FALSE
		#Get the TSSs from the annotation object
		chrs=as.character(seqnames(annot))
		strands=as.character(strand(annot))
		TSSs=start(annot)
		TSSs[strands=="-"]=end(annot)[strands=="-"]
		#Fetch the data for each thingy
		out=list()
		stop=length(annot)
		#w needs to be a whitelist, not a blacklist, or else things fail...
		w=c()
		for(i in 1:stop){
			if(strands[i]=="+"){
				tmp=getRegion(data,chrs[i],TSSs[i]-up,TSSs[i]+down)
				#Check for truncated genes
				if(length(tmp)!=up+down+1){
					out[[i]]="FAIL"
				}else{
					out[[i]]=tmp
					w=c(w,i)
				}
			}else{
				tmp=rev(getRegion(data,chrs[i],TSSs[i]-down,TSSs[i]+up))
				#Check for truncated genes
				if(length(tmp)!=up+down+1){
					out[[i]]="FAIL"
				}else{
					out[[i]]=tmp
					w=c(w,i)
				}
			}
			pp()
		}
		dat=out
		if(length(w)!=stop){
			message("Dropped ",stop-length(w)," genes due to inability to fetch data.")
			dat=dat[w]
		}
		xlabs=seq(-up,down,1)
	}else{
		dat=data
		row.ids=names(data)
		#Check that everything is nice and consistent (i.e., throw out bad genes)
		#First find the standard length
		st_len=max(sapply(dat,function(u){if("data"%in%names(u)){length(u$data)}else{0}}))
		#Identify bad genes
		w=which(sapply(dat,function(u){if("data"%in%names(u)){length(u$data)==st_len}else{FALSE}}))
		if(length(w)!=length(dat)){
			warning("Dropping ",length(dat)-length(w)," genes due to unavailable data.")
			dat=dat[w]
		}
		#Now we have to get the range requested, set it to be a constant subset.  We assume constant max/min/resolution
		drange=range(as.numeric(gsub("%","",dat[[1]]$samples))/dat[[1]]$glen)
		dres=dat[[1]]$resolution
		m=round(approx(seq(drange[1]*100,drange[2]*100,dres),1:length(dat[[1]]$samples),xout=seq(-up,down,dres))$y)
		xlabs=seq(drange[1]*100,drange[2]*100,dres)[m]
		dat=lapply(dat,function(u){u$data[m]})
	}
	#Now we have a list with each entry containing the oriented data for each gene, smooth if we have to
	#Make sure the smoothing is odd...
	smoothing=smoothing+1-smoothing%%2
	if(smoothing>1){
		message("Applying smoothing...")
		#smooth then sample
		#First smooth, but runmean crashes all the time soooo...
		tmp=lapply(dat,function(u){tryCatch((runmean(u,k=smoothing,endrule='constant')[seq(1,length(u),smoothing)]),error=function(code){NULL})})
		#If it's failed, try the more stringent thing
		fail=sapply(tmp,length)==0
		tmp[fail]=lapply(dat[fail],function(u){tryCatch((runmean(Rle(as.numeric(u)),k=smoothing,endrule='constant')[seq(1,length(u),smoothing)]),error=function(code){NULL})})
		#If there's any that still don't work, fuck em...
		fail=sapply(tmp,length)==0
		tmp[fail]=Rle(rep(0,3001))
		dat=tmp
		xlabs=xlabs[seq(1,length(xlabs),smoothing)]
		rm(tmp)
	}
	message("Converting to a matrix...")
	#Now we can convert it into a matrix
	y=matrix(unlist(lapply(dat,as.numeric)),nrow=length(dat),ncol=length(dat[[1]]),byrow=TRUE)
	#Cap it
	y_max=quantile(y,cap_quant)
	y[y>y_max]=y_max
	#Set the row.ids
	rownames(y)=row.ids[w]
	#Save it
	output$data=y

	#Perform the clustering
	if(cluster){
		message("Clustering...")
		cl=kmeans(y,ncluster,iter.max=100)
		output$cluster=cl
	}

	#What order should we plot stuff, up to a maximum of 5 things
	message("Sorting...")
	sorts=list()
	sorts[[1]]=sorts[[2]]=sorts[[3]]=sorts[[4]]=sorts[[5]]=rep(1,nrow(y))
	if(length(sort_by)!=0){
		for(i in length(sort_by):1){
			x=sort_by[i]
			rev=1
			if(substr(x,1,1)=="-"){
				rev=-1
				x=gsub("-",'',x)
			}
			if(x=="expression")
				sorts[[i]]=rev*exp_vector[w]
			if(x=="status")
				sorts[[i]]=rev*status_vector[w]
			if(x=="data")
				sorts[[i]]=rev*data_vector[w]
			if(x=="coverage")
				sorts[[i]]=rev*(rowSums(y)/ncol(y))
			if(x=="cluster")
				sorts[[i]]=rev*(cl$cluster)
		}
	}
	ordering=order(sorts[[1]],sorts[[2]],sorts[[3]],sorts[[4]],sorts[[5]])

	output$row.ids=row.ids[w]
	#If they are included, put them in the output list...
	if(!is.null(status_vector))
		output$status=status_vector[w]
	if(!is.null(exp_vector))
		output$expression=exp_vector[w]
	if(!is.null(data_vector))
		output$metadata=data_vector[w]
		
	#Perform the plot
	if(plot){
		message("Plotting...")
		#We define the layout...
		if(plot_key){
			ll=1:2
			ww=c(1.5,6)
		}else{
			ll=1
			ww=6
		}
		if(!is.null(status_vector)){
			ll=c(ll,length(ll)+1)
			ww=c(ww,1)
		}
		if(!is.null(exp_vector)){
			ll=c(ll,length(ll)+1)
			ww=c(ww,2)
		}
		if(!is.null(data_vector)){
			ll=c(ll,length(ll)+1)
			ww=c(ww,2)
		}
		layout(rbind(ll),widths=ww)
		par(mai=c(1.02,.5,.82,.05),oma=c(0,0,0,0))
		if(plot_key){
			#The signal intensity key
			image(rbind(1:64),col=cols,axes=FALSE,xlab="Enrichment\nLevel")
			#Put some axis on it
			axis(2,at=seq(0,1,length.out=10),labels=format(seq(0,y_max,length.out=10),digits=1))
		}
		#Plot the main thing
		image(xlabs,1:nrow(y),t(y[ordering,]),yaxt='n',xlab=xlab,ylab="",col=cols)
		if(cluster){
			#Add on break points for each cluster...
			br=cumsum(table(cl$cluster))
			for(b in br){lines(x=c(min(xlabs)-10,max(xlabs)+10),y=c(b,b),lwd=3)}
			#Add on a zero line
			lines(x=c(0,0),y=c(1,nrow(y)),lwd=3)
		}
		#Add on status
		if(!is.null(status_vector)){
			par(mai=c(1.02,.05,.82,.05))
			image(0,1:nrow(y),rbind(status_vector[w][ordering]),col=status_cols,axes=FALSE,xlab="Classification")
		}
		#Add on expression
		if(!is.null(exp_vector)){
			par(mai=c(1.02,.05,.82,.05))
			plot(exp_vector[w][ordering],1:nrow(y),yaxt='n',lwd=1,xlab="Expression",yaxs='i',pch='.',cex=4)
			qq=quantile(exp_vector[w],.75,na.rm=TRUE)
			#abline(v=qq,col=2,lty=2)
			#Add in something about the quantile average here
		}
		#Add on extra data
		if(!is.null(data_vector)){
			par(mai=c(1.02,.05,.82,.05))
			plot(data_vector[w][ordering],1:nrow(y),type='l',yaxt='n',lwd=1,yaxs='i')
			qq=quantile(data_vector[w],.75,na.rm=TRUE)
			abline(v=qq,col=2,lty=2)
		}
		#Finally, the title
		par(oma=c(0,0,2,0))
		mtext(title,line=0,outer=TRUE,cex=2)
	}
	invisible(output)
}

#Get rid of the factors in an efficient manner.  Have to do it on a case by case basis
unfactor=function(var){
	if (is.factor(var)){
		tmp=names(var)
		tmpopt=getOption("warn")
		options(warn=-1)
		out=as.numeric(levels(var))[as.integer(var)]
		options(warn=tmpopt)
		if(any(is.na(out)) & any(is.na(out)!=is.na(var))){
			out=as.character(levels(var))[as.integer(var)]
		}
		names(out)=tmp
	}else if(is.data.frame(var)){
		#Have to use a loop, since calling apply will return a matrix
		out=var
		for(i in 1:dim(var)[2]){
			out[,i]=unfactor(var[,i])
		}
	}else if(is.list(var)){
		#Mmmmm, recursion
		out=lapply(var,unfactor)
	}else{
		#The default option
		out=var
	}
	return(out)
}



#Prints the percentage done of a loop indexed by i given a total number of iterations stop
pp=function(total,count){
	if(missing(count)){count=evalq(i,envir=parent.frame())}
	if(missing(total)){total=evalq(stop,envir=parent.frame())}
	cat(round(100*(count/total)),"%   \r")
	if(count==total){cat("\n")}
}

#Function for classifying genes

classify_gene=function(gene,annot,cover,peak_window=200,min_gene_length=5000,Promoter_region=c(-3000,-100),TSS_region=c(-99,1000),Broad_region=c(1001,width(gene_info)),Promoter_ratio=1.25,TSS_ratio=1.25,Broad_ratio=.35){
	gene_info=annot[match(gene,values(annot)$gene_id)]
	#Is the gene too short?
	if(width(gene_info)<min_gene_length)
		return(list(classification="Unclassified",Promoter=c(NA,NA),TSS=c(NA,NA),Broad=c(NA,NA)))
	#We have two criteria to test, average coverage and peak height.  Start by calculating average coverage
	#What is the relevant region we need to get data for?
	data_region=c(min(Promoter_region,TSS_region,Broad_region)-peak_window,max(Promoter_region,TSS_region,Broad_region)+peak_window)
	forward=as.character(strand(gene_info))=="+"
	if(forward){
		fetch_region=start(gene_info)+data_region
	}else{
		fetch_region=rev(end(gene_info)-data_region)
	}
	#Get the actual data
	data=getRegion(cover,chr=as.character(seqnames(gene_info)),start=fetch_region[1],end=fetch_region[2])
	if(length(data)!=fetch_region[2]-fetch_region[1]+1)
		return(list(classification="Unclassified",Promoter=c(NA,NA),TSS=c(NA,NA),Broad=c(NA,NA)))
	if(!forward)
		data=rev(data)
	#Smooth it for peak determination
	peak_data=runmean(data,k=peak_window,endrule='constant')
	#Now we need to extract the three regions from it
	w=c(((Promoter_region[1]):(Promoter_region[2]))-data_region[1]+1)
	Promoter_data=data[w]
	Promoter_peak=peak_data[w]
	w=c(((TSS_region[1]):(TSS_region[2]))-data_region[1]+1)
	TSS_data=data[w]
	TSS_peak=peak_data[w]
	w=c(((Broad_region[1]):(Broad_region[2]))-data_region[1]+1)
	Broad_data=data[w]
	Broad_peak=peak_data[w]
	#Determine which region has the greatest coverage
	mean_cover=c(mean(Promoter_data),mean(TSS_data),mean(Broad_data))
	max_peak=c(max(Promoter_peak),max(TSS_peak),max(Broad_peak))
	max_cover=which(max(mean_cover)==mean_cover)
	classification="Unclassified"
	#The promoter region has the highest coverage
	if(length(max_cover)==1 & max_cover==1){
		if(max_peak[1]/max(max_peak[-1])>Promoter_ratio)
			classification="Promoter"
	}
	#The TSS region has the highest coverage
	if(length(max_cover)==1 & max_cover==2){
		if(max_peak[2]/max(max_peak[-2])>TSS_ratio)
			classification="TSS"
	}
	#The broad region has the highest coverage
	if(length(max_cover)==1 & max_cover==3){
		#How much of the signal is above the mean?  This removes broad genes with spikes
		if(sum(Broad_peak>mean(Broad_peak))/length(Broad_peak)>Broad_ratio)
			classification="Broad"
	}
	return(list(classification=classification,Promoter=c(max_peak[1],mean_cover[1]),TSS=c(max_peak[2],mean_cover[2]),Broad=c(max_peak[3],mean_cover[3])))
}

#############################################################################
#Description: For creating ASE plots.  Takes a coverage object, annotation object (range) and returns the averaged scaled signal from 5' to 3' of each annotation
#Notes: This is the simplified version where there is no constant bp flank and the smoothing is applied to the entire genome in a first step
#Author: Matthew Young
#Date Modified: 23/8/2010

getGenesScaledSimplified=function(cover,range,up=-100,down=200,res=1,stop=length(range)){
	samp=seq(up,down,res)
	n.bins=1+100/res
	#Pick a smoothing size that will use every bp of data even in the largest gene...
	max_gene=max(width(range))
	bw=max_gene/n.bins
	bw=bw-bw%%2+1
	#Now smooth the coverage object...
	smoothed=runmean(cover,bw,endrule="constant")
	bw.old=bw
	averaged=as.list(rep(NA,stop))
	chrs=as.character(seqnames(range))
	strands=as.character(strand(range))
	starts=as.integer(start(range))
	ends=as.integer(end(range))
	#starts<ends always, we want starts=5', ends=3'
	starts[strands=="-"]=ends[strands=='-']
	ends[strands=='-']=as.integer(start(range))[strands=='-']
	nom=as.character(values(range)[,1])
	#There's a weird bug in runmean, seems to be solved by rerunning a couple of times, this does that for us as much as is needed
	set=1:stop
	while(length(set)!=0){
		failed=rep(FALSE,stop)
		for(j in 1:length(set)){
			i=set[j]
			#We need to translate each samp point into a bp location in genomic space
			glen=abs(starts[i]-ends[i])+1
			#Where are we going to take as representative?  Remember, we're looking in a region ordered from upstream to downstream starting at index 1
			samples=floor(samp*glen*.01)
			offset=min(samples)-1
			samples=samples-offset
			i.up=abs(ceiling(up*glen*.01))
			i.down=abs(ceiling(down*glen*.01))+1
			if(strands[i]=="+"){
				reads=getRegion(smoothed,chr=chrs[i],start=starts[i]-i.up,end=starts[i]+i.down)
			}else{
				reads=rev(getRegion(smoothed,chr=chrs[i],start=starts[i]-i.down,end=starts[i]+i.up))
			}
			#If we didn't get all the bases we expected to, something went wrong, so skip this gene
			if(length(reads)!=i.up+i.down+1){
				averaged[[i]]="Data fetching error"
				warning("Was unable to get all the data for the genomic region around",starts[i])
				pp()
				next
			}
			#And now pick those points where we said we'd sample...
			data=seqselect(reads,start=samples,end=samples)
			#The % indicates which points are taken from the scaled central region
			coords=Rle(paste(samples+offset,"%",sep=''))
			#Push to output
			averaged[[i]]=list(data=data,samples=coords,smoothing=bw,resolution=res,glen=glen)
			pp()
		}
		set=which(failed)
	}
	names(averaged)=nom
	averaged
}

#############################################################################
#Description: Plots the average around the 5' end of an annotation with an up/down flank, or the ASE plot.
#Notes: If annot is given, then output is assumed to be a coverage object and data is fetched on the fly.  This only works for TSS plots, not ASE plots currently.
#range=TRUE, ASE plot, range=FALSE, TSS plot
#drop_truncated currently does nothing, the truncated genes are always dropped
#Author: Matthew Young
#Date Modified: 13/10/2010

plotAveraged = function(output,annot,upflank=10000,downflank=10000,lib.size=1e6,smoothing=1,blocksize=500,range=FALSE,add=FALSE,quartile=FALSE,invisible=FALSE,drop_truncated=TRUE,y_min=NULL,y_max=NULL,xlab="Bases from TSS",ylab="Average # Reads/bp/million",ylim=c(y_min,y_max),transform=function(u){u},...){
	if(xlab!=''){
		xlab="Bases from TSS"
		if(range)
			xlab="% gene length from TSS"
	}
	lib.size=lib.size*1e-6
	#If we have the annotation, we are fetching from data
	if(!missing(annot)){
		#Data fetching is only supported for TSS plots
		range=FALSE
		if(xlab!='')
			xlab="Bases from TSS"
		#Get the TSSs from the annotation object
		chrs=as.character(seqnames(annot))
		strands=as.character(strand(annot))
		TSSs=start(annot)
		TSSs[strands=="-"]=end(annot)[strands=="-"]
		#Fetch the data for each thingy
		out=rep(0,upflank+downflank+1)
		stop=length(annot)
		dropped=0
		for(i in 1:stop){
			if(strands[i]=="+"){
				tmp=transform(getRegion(output,chrs[i],TSSs[i]-upflank,TSSs[i]+downflank))
				#Check for truncated genes
				if(length(tmp)!=upflank+downflank+1){
					dropped=dropped+1
					next
				}
				if(smoothing>1){
					tmp=runmean(tmp,k=smoothing,endrule='constant')
				}
				if(quartile){
					out=cbind(out,tmp)
				}else{
					out=out+tmp
				}
			}else{
				tmp=rev(transform(getRegion(output,chrs[i],TSSs[i]-downflank,TSSs[i]+upflank)))
				#Check for truncated genes
				if(length(tmp)!=upflank+downflank+1){
					dropped=dropped+1
					next
				}
				if(smoothing>1){
					tmp=runmean(tmp,k=smoothing,endrule='constant')
				}
				if(quartile){
					out=cbind(out,tmp)
				}else{
					out=out+tmp
				}
			}
			pp()
		}
		if(quartile){
			quart=apply(out,1,boxplot,plot=FALSE)
			quart=sapply(quart,function(u){u$stats})
			out=apply(out,1,sum)
		}
		if(dropped>0){
			warning("Dropping ",dropped," genes due to incomplete data.")
		}
		if(smoothing>1){
			#out=runmean(out,k=smoothing,endrule='constant')
			if(quartile){
				quart=t(apply(quart,1,function(u){as.numeric(runmean(Rle(u),k=smoothing,endrule="constant"))}))
			}
		}
		#Make the xy thing
		xy=cbind(-upflank:downflank,as.numeric(out)/stop)
	}else{
		#We assume that the bad entries will have length 1, or have been filtered out beforehand
		w=which(sapply(output,length)==1)
		if(length(w)>0){
			warning("Dropping ",length(w)," genes due to unavailable data.")
			output=output[-w]
		}
		#Dropping the truncated genes
		tmp=table(sapply(output,function(u){length(u$data)}))
		master_len=as.numeric(names(tmp)[tmp==max(tmp)])
		w=which(sapply(output,function(u){length(u$data)})!=master_len)
		if(length(w)>0){
			warning("Dropping ",length(w)," genes due to incomplete data.")
			output=output[-w]
		}
		#We can now assume everything has the same length
		#If any of the samples things have a % in there, it has to be range
		if(length(grep("%",output[[1]]$samples))!=0)
			range=TRUE
		if(range){
			#First infer the start/end of the ranged region
			tmp=output[[1]]$samples
			#Is it the old type where there's no mixing between scaled and constant bp...
			all_scaled=FALSE
			if(length(grep("%",tmp))==0 | length(grep("%",tmp))==length(tmp))
				all_scaled=TRUE
			if(!all_scaled)
				tmp=tmp[grep("%",tmp)]
			master_range=range((as.numeric(gsub("%",'',tmp))-1-output[[1]]$glen)/output[[1]]$glen)
			#Now finesse it back to what it really is, because this involves a floor operation, we need to build a consensus
			start=table(sapply(output[1:min(1000,length(output))],function(u){
				tmp=u$samples[grep("%",u$samples)]
				if(length(tmp)==0)
					tmp=u$samples
				return(round(min(as.numeric(gsub("%",'',tmp)))/u$glen)*100)}))
			start=as.numeric(names(start)[start==max(start)])
			master_res=output[[1]]$resolution
			master_labels=seq(start,by=master_res,length.out=length(tmp))
			#Work out the number of bp up/down we are fetching
			if(!all_scaled){
				#The indicies tell us the sizes up/down
				tmp=grep("%",output[[1]]$samples)
				up.bp=min(tmp)-1
				down.bp=length(output[[1]]$samples)-length(tmp)-up.bp
			}
		}else{
			#We need to get the sampeled positions for each gene in our object
			master_range=range(output[[1]]$samples)
			master_labels=seq(master_range[1],master_range[2])
			#We want to pad everything to match the largest range given
			master=rep(0,length(master_labels))
		}
		if(quartile){
			if(!range | all_scaled){
				tt=sapply(output,function(u){as.numeric(u$data)})
			}else{
				tt=sapply(output,function(u){as.numeric(u$data[grep("%",u$samples)])})
			}
			m=apply(tt,1,mean)
			quart=apply(tt,1,boxplot,plot=FALSE)
			quart=sapply(quart,function(u){u$stats})
			if(smoothing>1){
				m=as.numeric(runmean(Rle(m),k=smoothing,endrule='constant'))
				quart=t(apply(quart,1,function(u){as.numeric(runmean(Rle(u),k=smoothing,endrule="constant"))}))
			}
			#We have to get flanks..
			if(range & !all_scaled){
				tt=sapply(output,function(u){as.numeric(u$data[grep("%",u$samples,invert=TRUE)])})
				flank.m=apply(tt,1,mean)
				flank.quart=apply(tt,1,boxplot,plot=FALSE)
				flank.quart=sapply(flank.quart,function(u){u$stats})
				if(smoothing>1){
					flank.m=as.numeric(runmean(Rle(flank.m),k=smoothing,endrule='constant'))
					flank.quart=t(apply(flank.quart,1,function(u){as.numeric(runmean(Rle(u),k=smoothing,endrule='constant'))}))
				}
			}
		}else{
			m=rep(0,length(master_labels))
			#Do this in blocks of blocksize, because otherwise pants are pooed, basically larger blocksize=more mem,more speed
			for(i in 1:ceiling(length(output)/blocksize)){
				if(!range | all_scaled){
					tmp=unlist(sapply(output[((i-1)*blocksize):min(length(output),i*blocksize)],function(u){as.numeric(u$data)}),FALSE,FALSE)
				}else{
					tmp=unlist(sapply(output[((i-1)*blocksize):min(length(output),i*blocksize)],function(u){as.numeric(u$data[grep("%",u$samples)])}),FALSE,FALSE)
				}
				m=apply(tmp,1,sum)+m
				pp(length(output),min(length(output),i*blocksize))
			}
			cat("\n")
			m=m/length(output)
			if(smoothing>1){
				m=as.numeric(runmean(Rle(m),k=smoothing,endrule='constant'))
			}
			#We have to get flanks...
			if(range & !all_scaled){
				flank.m=rep(0,up.bp+down.bp)
				#Do this in blocks of blocksize, because otherwise pants are pooed, basically larger blocksize=more mem,more speed
				for(i in 1:ceiling(length(output)/blocksize)){
					tmp=unlist(sapply(output[((i-1)*blocksize):min(length(output),i*blocksize)],function(u){as.numeric(u$data[grep("%",u$samples,invert=TRUE)])}),FALSE,FALSE)
					flank.m=apply(tmp,1,sum)+flank.m
					pp(length(output),min(length(output),i*blocksize))
				}
				cat("\n")
				flank.m=flank.m/length(output)
				if(smoothing>1){
					flank.m=as.numeric(runmean(Rle(flank.m),k=smoothing,endrule='constant'))
				}
			}
		}
		xy=cbind(master_labels,m)
		#Do we have to put flanks on these SOBs too?
		if(range & !all_scaled){
			y=c(head(flank.m,n=up.bp),m,tail(flank.m,n=down.bp))
			#Now we have to makeup some sensible labels...
			x=c(seq(-100+min(master_labels),min(master_labels),length.out=up.bp),master_labels,seq(max(master_labels),100+max(master_labels),length.out=down.bp))
			xy=cbind(x,y)
			#Might have to do the same business for the other quartiles...
			if(quartile){
				qquart=matrix(0,nrow=nrow(quart),ncol=ncol(quart)+ncol(flank.quart))
				for(i in 1:nrow(quart)){
					qquart[i,]=c(head(flank.quart[i,],n=up.bp),quart[i,],tail(flank.quart[i,],n=down.bp))
				}
				quart=qquart
			}
		}
	}	
	#Prepare the core data
	xy[,2]=xy[,2]/lib.size
	xy=xy[order(xy[,1]),]
	#Finally, do some plotting...
	if(!invisible){
		if(add){
			if(quartile){
				lines(xy[,1],xy[,2],lwd=3,...)
				lines(xy[,1],quart[2,]/lib.size,lwd=1,...)
				lines(xy[,1],quart[4,]/lib.size,lwd=1,...)
				lines(xy[,1],quart[1,]/lib.size,lwd=1,lty=3,...)
				lines(xy[,1],quart[5,]/lib.size,lwd=1,lty=3,...)
			}else{
				lines(xy[,1],xy[,2],...)
			}
		}else{
			#Construct y_min, y_max if needed
			if(length(ylim)==2){
				y_min=ylim[1]
				y_max=ylim[2]
			}
			#The user must specify all or none of the labels
			if(quartile){
				if(is.null(y_min))
					y_min=min(quart)/lib.size
				if(is.null(y_max))
					y_max=max(quart)/lib.size
				plot(xy[,1],xy[,2],type='l',xlab=xlab,ylab=ylab,ylim=c(y_min,y_max),lwd=3,...)
				lines(xy[,1],quart[2,]/lib.size,lwd=1,...)
				lines(xy[,1],quart[4,]/lib.size,lwd=1,...)
				lines(xy[,1],quart[1,]/lib.size,lwd=1,lty=3,...)
				lines(xy[,1],quart[5,]/lib.size,lwd=1,lty=3,...)
			}else{
				if(is.null(y_min))
					y_min=min(xy[,2])
				if(is.null(y_max))
					y_max=max(xy[,2])
				plot(xy[,1],xy[,2],type='l',xlab=xlab,ylab=ylab,ylim=c(y_min,y_max),...)
			}
			abline(v=0,col=2,lty=2)
			if(range){
				if(max(xy[,1])>=100){
					abline(v=100,col=2,lty=2)
				}
				if(!all_scaled){
					if(min(master_labels)!=0)
						abline(v=min(master_labels),col=2,lty=2)
					if(max(master_labels)!=100)
						abline(v=max(master_labels),col=2,lty=2)
				}
			}
		}
	}
	invisible(xy)
}


