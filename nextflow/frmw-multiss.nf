frmw_full = '~mario/framework/full/src/frmw_full'
frmw_quick = '~mario/framework/quick/src/frmw_quick'
frmw_kmers = '~mario/framework/kmers/src/frmw_kmers'
frmw_lcs = '~mario/framework/lcs/src/frmw_lcs'
frmw_upgma = '~mario/framework/upgma/src/frmw_upgma'
frmw_align = '~mario/framework/align/src/frmw_align'

process similarity_scores 
{
	input:
		tuple path(in), val(mode)
	output:
		path "${in.baseName}.in_tfa"
		path "${in.baseName}.${mode}.ss"

	script:
	print (in)
	if (mode == 'full')
		"""
		$frmw_full -ALIGN -INFILE=$in
		"""
	else if (mode == 'quick')
		"""
		$frmw_quick -ALIGN -INFILE=$in
		"""
	else if (mode == 'kmers')
		"""
		$frmw_kmers -ALIGN -INFILE=$in
		"""
	else
		"""
		$frmw_lcs -ALIGN -INFILE=$in
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
	infiles = Channel.fromPath(params.in + '*.in_tfa')
	ssmethods = Channel.of ('full', 'quick', 'kmers', 'lcs')
	params = infiles.combine (ssmethods)
	params.view ()
	(infile, ssfile) = similarity_scores (params)
	gtfile = upgma (infile, ssfile)
	alfile = alignment (infile, ssfile, gtfile)
}
