library(circlize)
library(ggplot2)

#-----------------------------img size------------------------------------
png_size <- function(matrix,fname) {
	if (nrow(matrix)>15)
		{ img = paste(fname,'circleplot.png',sep='_')
		p <- png(img,res=300, height=5000,width=5000)
		}
	else 
		{ img = paste(fname,'circleplot.png',sep='_')
		p <- png(img,res=300, height=3000,width=3000)
		}	
	return (p)
}

#-----------------------------circplot function---------------------------
circplot <- function(R,NR,gene) 
{ 
	#palette and png size definition
	pal_NR <- colorRampPalette(c("lightsalmon1","lightsalmon4"))
	pal_R <- colorRampPalette(c("palegreen1","palegreen4"))
	col.mat = c(pal_R(nrow(R)),pal_NR(nrow(NR)))
	fname <- paste(gene,'.png',sep='')
	mat <- rbind(R,NR)
	pal_grid <- colorRampPalette(c("gray0","gray90"))
	greys <- pal_grid(nrow(mat))
	names(greys) <- rownames(mat)
	png_size(mat,gene)
	#chorddiagram
	circos.par(canvas.xlim = c(-1.5, 1.5), canvas.ylim = c(-1.5, 1.5))
	chordDiagram(as.matrix(mat),col=col.mat, grid.col = greys, transparency = 0.5, annotationTrack = c("grid"),annotationTrackHeight = c(0.01, 0.01))
	circos.trackPlotRegion(track.index = 1, 
		panel.fun = function(x, y) {
	  		xlim = get.cell.meta.data("xlim")
	 		ylim = get.cell.meta.data("ylim")
	  		sector.name = get.cell.meta.data("sector.index")
	  		circos.text(mean(xlim), ylim[1]+2, sector.name, facing='clockwise', niceFacing = TRUE, adj =c(0, 0.5), col = "black",cex=1)}, 
		bg.border = NA) 
	circos.clear()
	legend("bottomleft", legend =c('Responder', 'Non Responder'), pch=16, cex=2,pt.cex=2, bty='n',col = c('palegreen', 'lightsalmon'))
	dev.off()
}

#----------------------------importing data-------------------------------

#alleles
R <- read.table('hla_R.tsv')
NR <- read.table('hla_NR.tsv')
A_R <- R[,1:2]
A_NR <- NR[,1:2]
B_R <- R[,3:4]
B_NR <- NR[,3:4]
C_R <- R[,5:6]
C_NR <- NR[,5:6]
circplot(A_R,A_NR,'A')
circplot(B_R,B_NR,'B')
circplot(C_R,C_NR,'C')

#supertypes
R <- read.table('hlaST_R.tsv')
NR <- read.table('hlaST_NR.tsv')
A_R <- R[,1:2]
A_NR <- NR[,1:2]
B_R <- R[,3:4]
B_NR <- NR[,3:4]
circplot(A_R,A_NR,'A_ST')
circplot(B_R,B_NR,'B_ST')
