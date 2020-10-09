#--------------------------importing libraries-----------------------------

library(ComplexHeatmap)
library(circlize)
library(ggplot2)

#-----------------------importing group data-------------------------------

info <- read.table("info.tsv", header=TRUE)
hla <- read.table("hla.tsv", header=TRUE)
gr_split <- split(info, info$group)
NR <- gr_split$NR$sample
R <- gr_split$R$sample

#-------------------------------palettes-----------------------------------

col_fun_r=colorRamp2(c(0,0.1,0.2,0.3,0.7,1),c("white","#EBE48B","#F2DD80","#DF8B5B","#CF7355","#A63945"))
col.hm=colorRamp2(c(0,1.5,3,4.5,6,8,10,15,20),c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"))
col_fun_ann=colorRamp2(c(0,0.3),c("white","turquoise4"))
pal <- colorRampPalette(c('ivory','orangered3'))
col_fun_h = colorRamp2(c(0,0.01,0.2,0.3,0.4,0.5,1),c(pal(7)))

#------------------------------setting img size----------------------------

png_size <- function(matrix,fname) {
	if (nrow(matrix)>30)
		{ img = paste(fname,'heatmap.png',sep='_')
		a <- png(img,res=300, height=7000,width=3000)
		}
	else 
		{ img = paste(fname,'heatmap.png',sep='_')
		a <- png(img,res=300, height=3000,width=2000)
		}
	return (a)
}
#--------------------------------HED---------------------------------------

#importing data
a<-read.table("HED.tsv",header=TRUE,row.name=1)
a=as.matrix(a)

#colorvector for HED groups
col <- character()
lab <- character()
fc <- function(data)
{for (i in data)
	{	
		{if (i %in% NR) 
			{n='azure'
			} 
		else 	
			{n='azure3'
			}
		} 
	col=c(col,n) 
	}
return(col)
}
fl <- function(data)
{for (i in data)
	{	
		{if (i %in% NR) 
			{g = 'NR'
			} 
		else 	
			{g = 'R'
			}
		} 
	lab=c(lab,g)
	}
return(lab)
}

#setting labels
col.ha=fc(rownames(a))
labels=fl(rownames(a))
names(col.ha) = labels

#annotations
ha = HeatmapAnnotation(distance = anno_boxplot(a,height = unit(3, "cm"), gp = gpar(fill = c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4"))))
ha2 = rowAnnotation(group=labels,col = list(group = col.ha),border=TRUE)

#heatmap
png_size(a,"HED")
Heatmap(a,
border = TRUE,
name = "Distance",
column_title = "HED", 
row_title = "Samples",
row_dend_width = unit(1.5, "cm"),
col = col.hm,
column_order = sort(colnames(a)),
top_annotation=ha,
left_annotation=ha2,
column_names_rot = 0,
column_names_gp = grid::gpar(fontsize = 10),
column_title_gp = gpar(fill = "lightyellow", col = "black", border ="black",fontsize=20),
na_col="black")

dev.off() 

#------------------------response frequencies----------------------------- 

#impoting data
c<-read.table("resp_a_freq.tsv",header=TRUE,row.name=1)
e<-read.table("resp_st_freq.tsv",header=TRUE,row.name=1)

#defining heatmap funtion
heatmap <- function(df) {
	f_diff<-as.matrix(df)[,4]
	df=as.matrix(df)[,2:3]
	#annotation	
	ha = rowAnnotation(f_diff = f_diff,
	gp = gpar(col = "black"),
	simple_anno_size = unit(1, "cm"), 
	col=list(f_diff=col_fun_ann))
	#heatmap
	a=Heatmap(df,
	border = TRUE,
	name = "Frequency", 
	column_title = "Allele frequency",
	col = col_fun_r,
	row_order = sort(rownames(df)),
	row_names_gp = grid::gpar(fontsize = 10),
	column_names_gp = grid::gpar(fontsize = 15),
	column_order = sort(colnames(df)),
	rect_gp = gpar(col = "white", lwd = 2),
	column_names_rot = 0,
	column_title_gp = gpar(fill = "lightyellow", col = "black", border = "black",fontsize=10))
	return(ha+a)
	}

png_size(c,'response_frequency')
heatmap(c)
dev.off()
png_size(e,'responseST_freq')
heatmap(e)
dev.off()


quit()
