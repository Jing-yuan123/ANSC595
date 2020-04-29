# ANSC595
pwd

#Load Modules

module load bioinfo mothur
module list

date +"%d %B %Y %H:%M:%S"

echo "Running mothur pipeline for my project "

#echo "We are testing our submission file"
#mothur "#make.contigs(file=gene.files, processors=12)"
#mothur "#summary.seqs(fasta=gene.trim.contigs.fasta);screen.seqs(fasta=gene.trim.contigs.fasta, group=gene.contigs.groups, maxambig=0, maxlength=275)"
#mothur "#get.current();summary.seqs(fasta=gene.trim.contigs.good.fasta);summary.seqs(fasta=current);summary.seqs()"
#mothur "#unique.seqs(fasta=gene.trim.contigs.good.fasta);count.seqs(name=gene.trim.contigs.good.names, group=gene.contigs.good.groups);summary.seqs(count=gene.trim.contigs.good.count_table)"
#mothur "#pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=12);rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta);summary.seqs(fasta=silva.v4.fasta)"
#mothur "#align.seqs(fasta=gene.trim.contigs.good.unique.fasta, reference=silva.v4.fasta,flip=t);summary.seqs(fasta=gene.trim.contigs.good.unique.align, count=gene.trim.contigs.good.count_table)"
#mothur "#screen.seqs(fasta=gene.trim.contigs.good.unique.align, count=gene.trim.contigs.good.count_table, summary=gene.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8);summary.seqs(fasta=current, count=current)"
#mothur "#filter.seqs(fasta=gene.trim.contigs.good.unique.good.align, vertical=T, trump=.);unique.seqs(fasta=gene.trim.contigs.good.unique.good.filter.fasta, count=gene.trim.contigs.good.good.count_table)"
#mothur "#pre.cluster(fasta=gene.trim.contigs.good.unique.good.filter.unique.fasta, count=gene.trim.contigs.good.unique.good.filter.count_table, diffs=2)"
#mothur "#chimera.vsearch(fasta=gene.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t);remove.seqs(fasta=gene.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos);summary.seqs(fasta=current, count=current)"
#mothur "#classify.seqs(fasta=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset16_022016.rdp.fasta, taxonomy=trainset16_022016.rdp.tax, cutoff=80);remove.lineage(fasta=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);summary.tax(taxonomy=current, count=current)"
#mothur "#dist.seqs(fasta=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)"
#mothur "#cluster(column=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"
#mothur "#make.shared(list=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)"
#mothur "#classify.otu(list=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy, label=1)"
#mothur "#rename.file(taxonomy=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)"
#mothur "#count.groups(shared=gene.opti_mcc.shared)"
#mothur "#sub.sample(shared=gene.opti_mcc.shared, size=11184)"
#mothur "#rarefaction.single(shared=gene.opti_mcc.shared, calc=sobs, freq=100)"
#mothur "#summary.single(shared=gene.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson-chao-shannon, subsample=T)"
#mothur "#dist.shared(shared=gene.opti_mcc.shared, calc=thetayc-jclass-braycurtis, subsample=t)"
#mothur "#pcoa(phylip=gene.opti_mcc.braycurtis.0.03.lt.ave.dist)"
#mothur "#nmds(phylip=gene.opti_mcc.braycurtis.0.03.lt.ave.dist,mindim=3, maxdim=3)"
#mothur "#nmds(phylip=gene.opti_mcc.braycurtis.0.03.lt.ave.dist)"
mothur "#amova(gene.opti_mcc.braycurtis.0.03.lt.ave.dist, design=ndesign.txt)"
mothur "#homova(gene.opti_mcc.braycurtis.0.03.lt.ave.dist, design=ndesign.txt)"
#mothur "#phylotype(taxonomy=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy)"
#mothur "#make.shared(list=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.tx.list, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=1)"
#mothur "#classify.otu(list=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.tx.list, count=gene.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy, label=1)"
#mothur "#dist.seqs(fasta=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, output=lt, processors=8)"
#mothur "#clearcut(phylip=gene.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)"
#mothur "#normalize.shared(shared=gene.opti_mcc.shared, method=totalgroup)"
#mothur "#tree.shared(shared=gene.final.an.unique_list.0.03.norm.shared, label=0.03, calc=braycurtis)"
#mothur "#metastats(shared=gene.opti_mcc.0.03.subsample.shared, design=mouse.time.design)"




echo " "
date +"%d %B %Y %H:%M:%S"
