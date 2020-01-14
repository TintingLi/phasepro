


plot_LCD_Motif_AAs <- function(motif.dt, lcd.dt, gene_symbol){
	#Usage:
	#	plot the AA distributions, the charges and the phase separatioin scores
	#of given gene.
	#parameters:
	#	dt: data.table contain the amino acid properties of a protein.
	#		col names of dt:
	#			amino_acid	FCR	NCPR	iupred2_score	achore2_score.
	#	gene_symbol: the gene used for analysis.
	#AA charges:	
	#positive: R(Arg)	K(Lys)	H(His)
	#negative: D(Asp)	E(Glu)
	#polar: S(Ser)	T(Thr)	N(Asn)	Q(Gln)
	#special: C(Cys)	G(Gly)	P(Pro)
	#Hydrophobic: A(Ala)	V(Val)	I(Ile)	L(Leu)	M(Met)	F(Phe)	Y(Tyr)	W(Trp)
	
	layout(matrix(c(1,1,1,2,3,4,5,6,7,8)), 10,1)
	par(mar =c(1,4,1,1), oma = c(4, 3, 2, 1))
	
	#AAs:----	
	plot(x =0, xlim=c(0,lcd.dt[,.N]),ylim=c(0,20), type = 'n', yaxt = 'n', 
		 ylab="amino acids", main = gene_symbol)
	#add y label:
	axis(2, at = seq(0.5, 20), tick = F, las =2, line = -0.7, 
		 labels = c('R','K','H', 
		 		   'D', 'E', 
		 		   'S','T', 'N','Q', 
		 		   'C', 'G', 'P',
		 		   'A', 'V', 'I', 'L','M','F', 'Y', 'W') %>% rev)
	
	#background color for the five classes of AAs.
	rect(0, 17, lcd.dt[, .N], 20, col = 1,  border = 1)
	rect(0, 15, lcd.dt[, .N], 17, col = 2,  border = 2)
	rect(0, 11, lcd.dt[, .N], 15, col = 3,  border = 3)
	rect(0, 8, lcd.dt[, .N], 11, col = 4,  border = 4)
	rect(0, 0, lcd.dt[, .N], 8, col = 5,  border = 5)
	
	
	aa_y_location <- data.table(
		amino_acid = c('R','K','H',  'D', 'E', 'S','T', 'N','Q', 'C', 
					   'G','P','A', 'V', 'I', 'L','M','F', 'Y', 'W'),
		location = seq(20,1))
	
	for (i in seq(1,lcd.dt[, .N])) {
		AA <- lcd.dt[i, amino_acid]
		j <- aa_y_location[amino_acid == AA, location]
		rect(i, j-1, i+1, j, col = par("fg"))
	}
	
	#Fraction of charged residue(FCR):----
	plot(1:lcd.dt[, .N], lcd.dt[, FCR], type = 'l', ylab = 'FCR', 
		 xaxt = 'n', xlab = '', lwd =2, col = 6)
	#Net charges per residue(NCPR): ----
	plot(lcd.dt[, NCPR], type = 'l', xaxt = 'n', ylab = 'NCPR',lwd =2, col = 7)
	abline(h = 0, col = "lightgray")	
	#Phase separation scores predicted by iupred2 software: ----
	plot(lcd.dt[, iupred2_score], type = 'l', ylab = 'phase separation',
		 xlab = '', col = 8, lwd =3, ylim = c(0,1))
	lines(lcd.dt[, achore2_score],type = 'l', col = 10, lwd =3)
	abline(h = 0.5, col = "lightgray")	
	
	# MOTIFS ------------------------------------------------------------------
	#pfam:----
	plot(x =0, xlim=c(0, motif.dt[1, AA_length ]), ylim=c(0,2), bty = 'n',
		 type = 'n', xaxt= 'n', yaxt = 'n', ylab="pfam")
	rect(0, 0.7, motif.dt[1, AA_length], 1.3, col = 1, border = 1)
	
	pfam.dt <- motif.dt[analysis_methods == 'Pfam']
	if (pfam.dt[, .N] == 0){
		mtext('no pfam motif found!')
	}else{
		for (i in seq(1, pfam.dt[, .N])) {
			x1 <- pfam.dt[i, start]
			x2  <- pfam.dt[i, end]
			rect(x1, 0.6, x2, 1.5, col = 2, border = 2)
			text((x1 + x2)/2, 1.85, labels = pfam.dt[i, motif])
		}
	}
	#LCD: ----
	lcd.dt <- motif.dt[analysis_methods == 'MobiDBLite']
	if (lcd.dt[, .N] > 0){
		for (i in seq(1, lcd.dt[, .N])) {
			x1 <- lcd.dt[i, start]
			x2  <- lcd.dt[i, end]
			rect(x1, 0.5, x2, 1.6, col = 3, border =3)
			text((x1 + x2)/2, 0.1, labels = lcd.dt[i, motif])
		}
	}
	#Coils: ----
	coil.dt <- motif.dt[analysis_methods == 'Coils']
	if (coil.dt[, .N] > 0){
		for (i in seq(1, coil.dt[, .N])) {
			x1 <- coil.dt[i, start]
			x2  <- coil.dt[i, end]
			rect(x1, 0.4, x2, 1.7, col = 4, border =4)
			text((x1 + x2)/2, 1, labels = coil.dt[i, motif])
		}
	}
	
	
	
	
	#smart:	----
	plot(x =0, xlim=c(0, motif.dt[1, AA_length ]), ylim=c(0,2), bty = 'n',
		 type = 'n', xaxt= 'n', yaxt = 'n', ylab="SMART")
	rect(0, 0.65, motif.dt[1, AA_length], 1.35, col = 1, border =1)
	
	SMART.dt <- motif.dt[analysis_methods == 'SMART']
	if (SMART.dt[, .N] > 0){
		for (i in seq(1, SMART.dt[, .N])) {
			x1 <- pfam.dt[i, start]
			x2  <- pfam.dt[i, end]
			rect(x1, 0.3, x2, 1.7, col = 5, border =5)
			text((x1 + x2)/2, 1.85, labels = SMART.dt[i, motif])
		}
	}
	
	#PRINTS ----
	plot(x =0, xlim=c(0, motif.dt[1, AA_length ]), ylim=c(0,2), bty = 'n',
		 type = 'n', xaxt= 'n', yaxt = 'n', ylab="PRINTS")
	rect(0, 0.65, motif.dt[1, AA_length], 1.35, col = 1, border =1)	
	
	PRINTS.dt <- motif.dt[analysis_methods == 'PRINTS']
	if (PRINTS.dt[, .N] > 0){
		for (i in seq(1, PRINTS.dt[, .N])) {
			x1 <- pfam.dt[i, start]
			x2  <- pfam.dt[i, end]
			rect(x1, 0.3, x2, 1.7, col = 6, border =6)
			text((x1 + x2)/2, 1.85, labels = PRINTS.dt[i, motif])
		}
	}
	
	#superfamily:	----
	plot(x =0, xlim=c(0, motif.dt[1, AA_length ]), ylim=c(0,2), bty = 'n',
		 type = 'n',yaxt = 'n', ylab="SUPERFAMILY")
	rect(0, 0.65, motif.dt[1, AA_length], 1.35, col = 1, border =1)	
	
	SUPERFAMILY.dt <- motif.dt[analysis_methods == 'SUPERFAMILY']
	if (SUPERFAMILY.dt[, .N] > 0){
		for (i in seq(1, SUPERFAMILY.dt[, .N])) {
			x1 <- pfam.dt[i, start]
			x2  <- pfam.dt[i, end]
			rect(x1, 0.3, x2, 1.7, col = 7, border =7)
			text((x1 + x2)/2, 1.85, labels = SUPERFAMILY.dt[i, motif])
		}
	}
	
	
	
	
	
	
	
	
	
	
	
}
