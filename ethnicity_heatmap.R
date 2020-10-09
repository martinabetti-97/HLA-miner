library(ComplexHeatmap)
library(circlize)
library(ggplot2)

#------------------------------setting img size----------------------------

png_size <- function(matrix,fname) {
	if (nrow(matrix)>30)
		{ img = paste(fname,'heatmap.png',sep='_')
		a <- png(img,res=300, height=7000,width=3000)
		}
	else 
		{ img = paste(fname,'heatmap.png',sep='_')
		a <- png(img,res=300, height=3000,width=1500)
		}
	return (a)
}

#-------------------------ethnicity frequencies----------------------------

#palette
col_fun_e=colorRamp2(c(0,0.1,0.2,0.3,0.4,0.5,0.7,1),c("white","#EBE48B","#F2DD80","#F6D277","#F6C56F","#DF8B5B","#CF7355","#A63945"))
	
#importing data
b<-read.table("ethnicity_freq.tsv",header=TRUE,row.name=1)
b=as.matrix(b)
#heatmap
png_size(b,'ethnicity_frequency')
Heatmap(b,
border = TRUE,
name = "Frequency", 
column_title = "Ethnic-specific frequencies",
col = col_fun_e,
column_dend_height = unit(2, "cm"),
row_order = sort(rownames(b)),
row_names_gp = grid::gpar(fontsize = 10), 
column_names_gp = grid::gpar(fontsize = 12), 
rect_gp = gpar(col = "black", lwd = 1),
column_title_gp = gpar(fill = "lightyellow", col = "black", border = "black",fontsize=15),column_names_rot = 0)
dev.off()

quit()
