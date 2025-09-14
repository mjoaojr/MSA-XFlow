frmw_full = '~/mnt/framework/full/src/frmw_full'
frmw_fullpar = '~/mnt/framework/fullpar/src/frmw_fullpar'
frmw_quick = '~/mnt/framework/quick/src/frmw_quick'
frmw_kmers = '~/mnt/framework/kmers/src/frmw_kmers'
frmw_lcs = '~/mnt/framework/lcs/src/frmw_lcs'
frmw_proba = '~/mnt/framework/proba/src/frmw_proba'
frmw_cons = '~/mnt/framework/consistency/src/frmw_cons'
frmw_upgma = '~/mnt/framework/upgma/src/frmw_upgma'
frmw_slmin = '~/mnt/framework/slmin/src/frmw_slmin'
frmw_slmax = '~/mnt/framework/slmax/src/frmw_slmax'
frmw_nj = '~/mnt/framework/nj/src/frmw_nj'
frmw_align = '~/mnt/framework/align/src/frmw_align'
mclustalw = '~/mnt/clustalw/mclustalw/src/clustalw2'

process full 
{
	maxcpus = Runtime.runtime.availableProcessors()
	cpus = maxcpus - 3
	input:
		path(in)
	output:
		tuple path ("${in.baseName}.in_tfa"), path ("${in.baseName}.${mode}.ss"), val (mode)

	script:
		mode = 'fullpar'
		"""
		$frmw_fullpar -ALIGN -INFILE=$in -NTHR=$cpus
		"""
}

process similarity_scores 
{
	input:
		tuple path(in), val(mode)
	output:
		tuple path ("${in.baseName}.in_tfa"), path ("${in.baseName}.${mode}.ss"), val (mode)

	script:
//	print ('SS: ' + in + ' ' + mode)
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
	else if (mode == 'proba')
		"""
		$frmw_proba -ALIGN -INFILE=$in
		"""
	else if (mode == 'cons')
		"""
		$frmw_cons -ALIGN -INFILE=$in
		"""
	else
		"""
		$frmw_lcs -ALIGN -INFILE=$in
		"""
}

process guide_tree 
{
	input:
		tuple path (infile), path (ss), val (ssmode), val (mode)
	output:
		tuple path (infile), path (ss), path ("${infile.baseName}.dnd"), val (ssmode), val (mode)

	script:
//	print ('GT: ' + infile + ' ' + ss + ' ' + ssmode + ' ' + mode)
	if (mode == 'upgma')
		"""
		$frmw_upgma -ALIGN -INFILE=$infile -SSFILE=$ss
		"""
	else if (mode == 'slmin')
		"""
		$frmw_slmin -ALIGN -INFILE=$infile -SSFILE=$ss
		"""
	else if (mode == 'slmax')
		"""
		$frmw_slmax -ALIGN -INFILE=$infile -SSFILE=$ss
		"""
	else
		"""
		$frmw_nj -ALIGN -INFILE=$infile -SSFILE=$ss
		"""
}

process alignment 
{
	input:
		tuple path (infile), path (ss), path (tree), val (ssmode), val (gtmode)
	output:
		tuple path (infile), path (ss), path (tree), val (ssmode), val (gtmode), path ("${infile.baseName}.msf")

	script:
//	print ('Alignment: ' + infile + ' ' + ss + ' ' + tree + ' ' + ssmode + ' ' + gtmode)
	gtopt = gtmode.toUpperCase()
	"""
	$frmw_align -ALIGN -CLUSTERING=$gtopt -INFILE=$infile -SSFILE=$ss -USETREE=$tree -MAXDIV=0 -OUTPUT=GCG
	"""
}

process check_align 
{
	input:
		tuple path (infile), path (ss), path (tree), val (ssmode), val (gtmode), path (alfile)
	output:
		path alfile

	script:
	outfile = infile.baseName + '.' + ssmode + '.' + gtmode + '.msf'
	ssopt = ssmode.toUpperCase()
	gtopt = gtmode.toUpperCase()
//	print ('Check: ' + infile + ' ' + ss + ' ' + tree + ' ' + ssmode + ' ' + gtmode + ' ' + outfile)
	if (ssmode == 'full')
		"""
		$mclustalw -ALIGN -CLUSTERING=$gtopt -INFILE=$infile -MAXDIV=0 -OUTPUT=GCG -OUTFILE=$outfile
		echo diff $outfile $alfile
		"""
	else if (ssmode == 'quick')
		"""
		$mclustalw -ALIGN -QUICKTREE -CLUSTERING=$gtopt -INFILE=$infile -MAXDIV=0 -OUTPUT=GCG -OUTFILE=$outfile
		echo diff $outfile $alfile
		"""
	else
		"""
		$mclustalw -ALIGN -$ssopt -CLUSTERING=$gtopt -INFILE=$infile -MAXDIV=0 -OUTPUT=GCG -OUTFILE=$outfile
		echo diff $outfile $alfile
		"""
}

workflow 
{
	infiles = Channel.fromPath(params.in + '*.in_tfa')
	infiles_full = Channel.fromPath(params.in + '*.in_tfa')
//	ssmethods = Channel.of ('full', 'quick', 'kmers', 'lcs', 'proba', 'cons')
	ssmethods = Channel.of ('quick', 'kmers', 'lcs')
	params = infiles.combine (ssmethods)
//	params.view ()
	fullss = full (infiles_full)
	ssfiles = similarity_scores (params)
	ssfiles.mix (fullss)
	gtmethods = Channel.of ('upgma', 'slmin', 'slmax', 'nj')
	params2 = ssfiles.combine (gtmethods)
//	params2.view ()
	gtfile = guide_tree (params2)
	alfile = alignment (gtfile)
//	alfinal = check_align (alfile)
}
