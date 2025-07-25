## $Id: Makefile.am 1421 2006-11-01 01:15:40Z cegrant $
##
## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = no-dependencies

SUBDIRS = motif.crp0 motif.Puf3p

DATAFILES = \
	At.s \
	DNaseHS_gm12878.wig \
	GCN4.meme.txt \
	GCN4_YPD-genomic.fasta \
	GCN4_YPD-genomic.psp \
	GCN4_YPD.fasta \
	GCN4_YPD.psp \
	GCN4_YPD.wig \
	GO2Gene.map.csv \
	INO_up800.s \
	JASPAR_CORE_2014_vertebrates.meme \
	Jaspar-subset.meme \
	Jaspar-subset2.meme \
	Klf1-200-100-shuffled.fa \
	Klf1-200-100-shuffled.py2.fa \
	Klf1-200-100-shuffled.py3.fa \
	Klf1-200-100.combined.fa \
	Klf1-200-100.combined.nofasta.fa \
	Klf1-200-100.fa \
	Klf1-200-100.no_fasta.fa \
	Klf1-200-100_ucsc_galaxy.fa \
	Klf1-200.fa \
	Klf1.100.fa \
	Klf1.fa \
	MCM1.fasta \
	MCM1.meme.html \
	Puf3p-20.s \
	Puf3p-200.s \
	Puf3p-dna_encoded.s \
	Puf3p.s \
	a.bed \
	adh.s \
	all-A.fasta \
	alphabet_dna.txt \
	alphabet_rna.txt \
	ama.nozscoring.xml \
	ama.zscoring.xml \
	b.bed \
	chr22.fasta \
	chrX.fa \
	crp0-test.fasta \
	crp0.fasta \
	crp0.meme.html \
	crp0.meme.txt \
	crp0.meme.xml \
	crp0.s \
	dna_with_covalent.alph \
	dna_with_covalent.bg \
	dna_with_covalent.fa \
	dna_with_covalent.meme \
	dreme-Klf1.xml \
	dros.fasta \
	dros.oops.html \
	ensembl_Plasmodium_falciparum_38.200.fa \
	extended_dna.alph \
	farntrans5.s \
	gomo_motifs.meme.txt \
	gomo_seqs.bg \
	gomo_seqs.fasta \
	gomo_seqs.norc.bg \
	gomo_seqs_red.fasta \
	gomoscore-motif23.xml \
	gomoscore.xml \
	lex0.s \
	limits.fasta \
	lipo-test.fasta \
	lipo.fasta \
	lipo.meme.html \
	lipo.meme.old.html \
	lipocalin.s \
	mcast_motifs.meme \
	mcast_seqs.fa \
	meme.lex0.zoops.txt \
	mini-drosoph.s \
	motif-scan.gff \
	motif.Puf3p \
	motif.crp0 \
	non-uniform.bg \
	noq.s \
	nr.freq \
	nrmix.100bp.fa \
	nt.freq \
	null.txt \
	observed.txt \
	one-peak-dna.fasta \
	one-peak-protein.fasta \
	primary.meme \
	prior.dist.txt \
	psp_test.psp \
	psp_test.s \
	random.fasta \
	sample.meme \
	secondary.meme \
	simple.aln \
	simple.meme \
	small-comet-mod.txt \
	small-msgf-mod.txt \
	small-tide-mod-1-out-of-2.txt \
	small-tide-mod-2-out-of-2.txt \
	small-tide-mod-phospho.fasta \
	small-tide-mod-phospho.txt \
	small-tide-mod.txt \
	small-yeast-prealigned-no-unknowns.txt \
	small-yeast.fasta \
	spiked.aln \
	spiked.fasta \
	test.aln \
	test_tree.txt \
	text-in-weights.ptm.txt \
	tiny_genome.fa \
	tiny_genome.fa.fai \
	tiny_mouse.fasta \
	uniform.txt \
	uprot.freq \
	yeast.nc.1.freq \
	yeast.nc.3.freq \
	yeast.nc.6.freq \
	yeast.tree

if HAVE_PYTHON_3
Klf1-200-100-shuffled.fa: Klf1-200-100-shuffled.py3.fa
	cp Klf1-200-100-shuffled.py3.fa $(top_builddir)/tests/common/Klf1-200-100-shuffled.fa
else
Klf1-200-100-shuffled.fa: Klf1-200-100-shuffled.py2.fa
	cp Klf1-200-100-shuffled.py2.fa $(top_builddir)/tests/common/Klf1-200-100-shuffled.fa
endif

check_SCRIPTS = Klf1-200-100-shuffled.fa

EXTRA_DIST = $(DATAFILES)
