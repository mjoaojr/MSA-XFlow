frmw_full = '~/mnt/framework/full/src/frmw_full'
frmw_quick = '~/mnt/framework/quick/src/frmw_quick'
frmw_kmers = '~/mnt/framework/kmers/src/frmw_kmers'
frmw_lcs = '~/mnt/framework/lcs/src/frmw_lcs'
frmw_proba = '~/mnt/framework/proba/src/frmw_proba'
frmw_cons = '~/mnt/framework/consistency/src/frmw_cons'
frmw_conspar = '~/mnt/framework/conspar/src/frmw_conspar'
frmw_upgma = '~/mnt/framework/upgma/src/frmw_upgma'
frmw_slmin = '~/mnt/framework/slmin/src/frmw_slmin'
frmw_slmax = '~/mnt/framework/slmax/src/frmw_slmax'
frmw_nj = '~/mnt/framework/nj/src/frmw_nj'
frmw_align = '~/mnt/framework/align/src/frmw_align'
frmw_aligncons = '~/mnt/framework/aligncons/src/frmw_aligncons'
mclustalw = '~/mnt/clustalw/mclustalw/src/clustalw2'

process full 
{
	input:
		path(in)
	output:
		tuple val (in.baseName), path ("${in.baseName}.tfa"), path ("${in.baseName}.${mode}.ss"), val (mode)

	script:
		mode = 'full'
		"""
		$frmw_full -ALIGN -INFILE=$in
		"""
}

process consistency 
{
	maxForks 1
	cpus = Runtime.runtime.availableProcessors() - params.n
//	cpus = 2 - params.n
	input:
		path(in)
	output:
		tuple val (in.baseName), path ("${in.baseName}.tfa"), path ("${in.baseName}.${mode}.ss"), val (mode), path ("${in.baseName}.cl")

	script:
//		print (task.cpus)
		mode = 'proba'
		"""
		$frmw_conspar -ALIGN -INFILE=$in -NTHR=${task.cpus}
		"""
}

process similarity_scores 
{
	maxForks params.n
	input:
		tuple path(in), val(mode)
	output:
		tuple val (in.baseName), path ("${in.baseName}.tfa"), path ("${in.baseName}.${mode}.ss"), val (mode)

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
	else
		"""
		$frmw_lcs -ALIGN -INFILE=$in
		"""
}

process guide_tree 
{
	input:
		tuple val (infilename), path (infile), path (ss), val (ssmode), val (mode)
	output:
		tuple val (infilename), path (infile), path (ss), path ("${infile.baseName}.dnd"), val (ssmode), val (mode)

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
		tuple val (filename), path (infile), path (ss), path (tree), val (ssmode), val (gtmode)
	output:
		tuple val (filename), path (infile), path (ss), path (tree), val (ssmode), val (gtmode), val (clmode), path ("${infile.baseName}.msf")

	script:
//	print ('Alignment: ' + infile + ' ' + ss + ' ' + tree + ' ' + ssmode + ' ' + gtmode)
	gtopt = gtmode.toUpperCase()
	clmode = 0
	"""
	$frmw_align -ALIGN -CLUSTERING=$gtopt -INFILE=$infile -SSFILE=$ss -USETREE=$tree -MAXDIV=0 -OUTPUT=GCG
	"""
}

process alignment_cons 
{
	input:
		tuple val (filename), path (infile), path (ss), path (tree), val (ssmode), val (gtmode), path (constraints)
	output:
		tuple val (filename), path (infile), path (ss), path (tree), val (ssmode), val (gtmode), val (clmode), path ("${infile.baseName}.msf")

	script:
//	print ('Alignment: ' + infile + ' ' + ss + ' ' + tree + ' ' + ssmode + ' ' + gtmode)
	gtopt = gtmode.toUpperCase()
	clmode = 1
	"""
	$frmw_aligncons -ALIGN -CLUSTERING=$gtopt -INFILE=$infile -SSFILE=$ss -USETREE=$tree -CLFILE=$constraints -MAXDIV=0 -OUTPUT=GCG
	"""
}

process copia
{
	input:
		tuple val (filename), path (infile), path (ss), path (tree), val (ssmode), val (gtmode), val (clmode), path (msf)
		val out
	output:
		path msf

	script:
//	print ('Alignment: ' + infile + ' ' + ss + ' ' + tree + ' ' + ssmode + ' ' + gtmode)
	gtopt = gtmode.toUpperCase()
	if (clmode == 1)
		outfile = out + filename + '.' + ssmode + '.' + gtmode + '.cons.msf'
	else
		outfile = out + filename + '.' + ssmode + '.' + gtmode + '.msf'
//	print (msf + ' to ' + outfile)
	xml = filename + '.xml'
	"""
	/usr/bin/cp $msf $outfile
        /home/ubuntu/mnt/BAliBASE_R9/src/bali_score ~/mnt/BAliBASE/RV100/$xml $outfile > eca
        echo $filename > scores
        echo $ssmode >> scores
        echo $gtmode >> scores
        echo $clmode >> scores
        cat eca | grep SP | cut -d' ' -f3 >> scores
        cat eca | grep TC | cut -d' ' -f3 >> scores

	python3 ~/mnt/savescores3.py scores BAliBASE RV100
	"""
//        /home/ubuntu/mnt/bali_score_src/bali_score ~/mnt/BAliBASE/RV100/$xml $outfile > eca
//        cat eca | grep CS | cut -d' ' -f3 >> scores
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
	infiles = Channel.fromPath(params.in + '*.tfa')
//	infiles_full = infiles
	infiles_cons = infiles
//	infiles_full = Channel.fromPath(params.in + '*.tfa')
//	ssmethods = Channel.of ('full', 'quick', 'kmers', 'lcs', 'proba', 'cons')
	ssmethods = Channel.of ('full', 'quick', 'kmers', 'lcs')
	parametros = infiles.combine (ssmethods)
//	params.view ()
//	infiles_full.view()
//	fullss = full (infiles_full)
//	fullss.view ()
	constraints = consistency (infiles_cons)
	ssfiles = similarity_scores (parametros)
//	ssfiles = ssfiles.mix (fullss)
	consssfiles = constraints.map { [it[0], it[1], it[2], it[3]] }
	ssfiles = ssfiles.mix (consssfiles)
//	ssfiles.view ()
	gtmethods = Channel.of ('upgma', 'slmin', 'slmax', 'nj')
	params2 = ssfiles.combine (gtmethods)
//	params2.view ()
	gtfile = guide_tree (params2)
	constraints = constraints.map { [it[0], it[4]] }
	gtconsfiles = gtfile.combine (constraints, by: 0)
	alfiles = alignment (gtfile)
//	gtfile.view()
//	alfinal = check_align (alfile)
//	gtconsfiles.view ()
	alconsfiles = alignment_cons (gtconsfiles)
	msffiles = alfiles.mix (alconsfiles)
	finalfiles = copia (msffiles, params.out)
}
