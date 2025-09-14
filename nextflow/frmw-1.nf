frmw_full = '~mario/framework/full/src/frmw_full'
frmw_upgma = '~mario/framework/upgma/src/frmw_upgma'
frmw_align = '~mario/framework/align/src/frmw_align'

process full 
{
	input:
		path in
	output:
		path "${in.baseName}.in_tfa"
		path "${in.baseName}.full.ss"

	script:
	"""
	$frmw_full -ALIGN -INFILE=$in
	"""
}

process upgma 
{
	input:
		path infile
		path ss
	output:
		path "${infile.baseName}.dnd"

	script:
	"""
	$frmw_upgma -ALIGN -INFILE=$infile -SSFILE=$ss
	"""
}

process alignment 
{
	input:
		path infile
		path ss
		path tree
	output:
		path "${infile.baseName}.msf"

	script:
	"""
	$frmw_align -ALIGN -INFILE=$infile -SSFILE=$ss -USETREE=$tree -MAXDIV=0 -OUTPUT=GCG
	"""
}

workflow 
{
	(infile, ssfile) = full (params.in)
	gtfile = upgma (infile, ssfile)
	alfile = alignment (infile, ssfile, gtfile)
}
