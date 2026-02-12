About the dataset…
================
Margaux Lefebvre
2026-02-12

# Sources

The samples are drawn from several literature sources:

- [Ibrahim *et al.* (2024)](https://doi.org/10.1038/s41467-024-55102-3)
  (n=222 *P. malariae*)
- [Plenderleith *et al*
  (2022)](https://doi.org/10.1038/s41467-022-29306-4) (n=1 *P. malariae*
  and n=1 *P. malariae-like*)
- [Ibrahim *et al* (2020)](10.1038/s41598-020-67568-4) (n=17 *P.
  malariae*)
- [Rutledge *et al* (2017)](https://doi.org/10.1038/nature21038) (n=4
  *P. malariae* and n=2 *P. malariae-like*)
- [Ansari *et al* (2016)](https://doi.org/10.1016/j.ijpara.2016.05.009)
  (n=1 *P. malariae*)

We also added 59 *P. malariae* and 20 *P. brasilianum* newly sequenced
isolates.

# Unfiltered BAMs and VCF with all the samples

## Mapping

Version: cutadapt v1.18, bwa-mem v0.7.17, samtools v1.15.1.

``` bash
echo "--> Start processing: $downId" #downId is the name of the sample

#Rename the files in order to remain consistent with the rest of the script
mv ${downId}_1.fastq.gz ${downId}.R1.fastq.gz
mv ${downId}_2.fastq.gz ${downId}.R2.fastq.gz

zcat ${downId}.R1.fastq.gz | sed -E 's/^((@|\+)'$downId'\.[^.]+)\.(1|2)/\1/' | bgzip > $downId.raw.1.fastq.gz
zcat ${downId}.R2.fastq.gz | sed -E 's/^((@|\+)'$downId'\.[^.]+)\.(1|2)/\1/' | bgzip > $downId.raw.2.fastq.gz

# Remove adapters and preprocessed to eliminate low-quality reads. Reads shorter than 50 bp containing “N” were discarded.
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -q 30 -m 25 --max-n 0 -o $path_fasta/$downId.R1.fastq.gz -p $path_fasta/$downId.R2.fastq.gz $path_fasta/$downId.raw.1.fastq.gz $path_fasta/$downId.raw.2.fastq.gz

#Sequenced reads were aligned to the P. malariae reference genome PmUG01
bwa mem -t 1 GCF_900090045.1_PmUG01_genomic.fna $downId.R1.fastq.gz $downId.R2.fastq.gz | samtools view -F 4 -b - | samtools sort - -o $downId.mapPM.sort.bam
samtools index $downId.mapPM.sort.bam
```

## Calling

Version: GATK v3.8.0, samtools v1.15.1, Picard tools v2.23.5.

``` bash
# Create Sequence Dictionary
picard CreateSequenceDictionary R=GCF_900090045.1_PmUG01_genomic.fna O=GCF_900090045.1_PmUG01_genomic.dict

# Mark the duplicated reads
picard MarkDuplicates -Xmx150g INPUT=$downId.mapPM.sort.bam OUTPUT=$downId.map.sort.dedup.bam METRICS_FILE=$downId.metrics.txt
samtools index $downId.map.sort.dedup.bam

#Add or Replace Read Groups
picard AddOrReplaceReadGroups -Xmx150g I=$downId.map.sort.dedup.bam O=$downId.map.sort.dedup.rg.bam \
LB=LIB-$downId PL=ILLUMINA PU=H0164ALXX140820:2:1101 SM=$downId

rm $downId.map.sort.dedup.bam $downId.map.sort.dedup.bam.bai
samtools index $downId.map.sort.dedup.rg.bam

# SplitNCigarReads
GenomeAnalysisTK -Xmx150g -T SplitNCigarReads -R GCF_900090045.1_PmUG01_genomic.fna -I $downId.map.sort.dedup.rg.bam \
-o $downId.map.sort.dedup.rg.rq.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

rm $bam_dir/$downId.map.sort.dedup.rg.bam $bam_dir/$downId.map.sort.dedup.rg.bam.bai
samtools index $bam_dir/$downId.map.sort.dedup.rg.rq.bam

# Local realignment around indels
GenomeAnalysisTK -Xmx150g -T RealignerTargetCreator -R GCF_900090045.1_PmUG01_genomic.fna \
-I $downId.map.sort.dedup.rg.rq.bam -o  $downId.realignertargetcreator.intervals

GenomeAnalysisTK -Xmx150g -T IndelRealigner -R GCF_900090045.1_PmUG01_genomic.fna -targetIntervals $downId.realignertargetcreator.intervals \
-I  $downId.map.sort.dedup.rg.rq.bam -o $downId.map.sort.dedup.indelrealigner.bam
mv $downId.map.sort.dedup.indelrealigner.bam $downId.bwa.gatk.sort.bam

samtools index $downId.bwa.gatk.sort.bam

# Calling with HaplotypeCaller
GenomeAnalysisTK -Xmx150g -T HaplotypeCaller -R GCF_900090045.1_PmUG01_genomic.fna \
-I $downId.bwa.gatk.sort.rescaled.bam --genotyping_mode DISCOVERY -stand_call_conf 10 \
-o $downId.ploidy2.raw_variants.snp.indel.g.vcf -ERC GVCF --sample_ploidy 2
```

## Combining all the samples and keep only the nuclear genome

Version: GATK v3.8.0, bcftools v1.16.

``` bash
#Merge all the samples
vcf_list=$(cat vcf_filtered.ploidy2.lst | awk '{print "--variant " $0}' | tr '\n' ' ')

GenomeAnalysisTK -Xmx50g -T CombineGVCFs \
  -R GCF_900090045.1_PmUG01_genomic.fna \
  $vcf_list \
  -o malariae_brasi_unfiltered.ploidy2.vcf

bgzip malariae_brasi_unfiltered.ploidy2.vcf
tabix malariae_brasi_unfiltered.ploidy2.vcf.gz

#Keep only the variants
GenomeAnalysisTK -Xmx50g -T GenotypeGVCFs -R GCF_900090045.1_PmUG01_genomic.fna \
-V malariae_brasi_total.ploidy2.vcf.gz -o  malariae_brasi.snps.ploidy2.vcf 
bgzip  malariae_brasi.snps.ploidy2.vcf
tabix  malariae_brasi.snps.ploidy2.vcf.gz

#Keep nuclear core_genome
bcftools view -R core_malariae.bcf.txt -o malariae_brasi.snps.core.ploidy2.vcf.gz malariae_brasi.snps.ploidy2.vcf.gz
tabix malariae_brasi.snps.core.ploidy2.vcf.gz 
```

# Filters on the dataset

## Quality and depth filtering

Version: bcftools v1.16, vcftools v0.1.16.

``` bash
VCF=malariae_brasi.snps.core.ploidy2.vcf.gz
OUT=./QC_VCF/malariae_brasi.info.ploidy2

vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
vcftools --gzvcf $VCF --depth --out $OUT
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
vcftools --gzvcf $VCF --site-quality --out $OUT
vcftools --gzvcf $VCF --missing-indv --out $OUT 
vcftools --gzvcf $VCF --missing-site --out $OUT
```

Filtering info:

- The minimum variant quality (Phred score) is 10, so we have to filter
  to 20.
- For the variant mean depth, we will set the minimum at 7X. The maximum
  will be set at 84X (= mean depth + twice the standard-deviation).
- For the individual mean depth, we will set the minimum at 7X and the
  maximum at 182X (= mean depth + twice the standard-deviation).
- We remove all the individuals with more than 75% of missing data: it
  removes 69 samples.
- Remove all the SNPs with more than 50% of missing data.
- To avoid sequencing error, we usually put the MAF at 1/number of
  samples : 1/327~0.003.

Applying filters to VCF

``` bash
VCF_IN=malariae_brasi.snps.core.ploidy2.vcf.gz 
VCF_OUT=malariae_brasi.filtsnps.core.ploidy2.vcf.gz

# set filters
MAF=0.003
MISS=0.5
QUAL=30
MIN_DEPTH_SNP=10
MAX_DEPTH_SNP=84
MIN_DEPTH=10
MAX_DEPTH=182

vcftools --gzvcf $VCF_IN --remove  missing_rm.samples --min-alleles 2 --max-alleles 2 \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH_SNP --max-meanDP $MAX_DEPTH_SNP \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > $VCF_OUT
```

What have we done here?

- `--remove remove_miss.txt` - remove all all the individuals with more
  than 50% of missing data
- `--min-alleles 2 --max-alleles 2` - keep only the bi-allelic SNPs
- `--remove-indels` - remove all indels (SNPs only)
- `--maf`- set minor allele frequency - here 1/number of samples
- `--max-missing` - set minimum missing data. A little counter
  intuitive - 0 is totally missing, 1 is none missing. Here 0.5 means we
  will tolerate 50% missing data.
- `--minQ` - this is just the minimum quality score required for a site
  to pass our filtering threshold. Here we set it to 20.
- `--min-meanDP` - the minimum mean depth for a site.
- `--max-meanDP` - the maximum mean depth for a site.
- `--minDP` - the minimum depth allowed for a genotype - any individual
  failing this threshold is marked as having a missing genotype.
- `--maxDP` - the maximum depth allowed for a genotype - any individual
  failing this threshold is marked as having a missing genotype.

## Remove multi-clonal infections (*F<sub>WS</sub>*)

> The F<sub>WS</sub> metric estimates the heterozygosity of parasites
> (HW) within an individual relative to the heterozygosity within a
> parasite population (HS) using the read count of alleles.
> F<sub>WS</sub> metric calculation for each sample was performed using
> the following equation: F<sub>WS</sub>=1− HW/HS where HW refers to the
> allele frequency of each unique allele found at specific loci of the
> parasite sequences within the individual, and HS refers to the
> corresponding allele frequencies of those unique alleles within the
> population. F<sub>WS</sub> ranges from 0 to 1; a low F<sub>WS</sub>
> value indicates low inbreeding rates within the parasite population
> and thus high within-host diversity relative to the population. An
> F<sub>WS</sub> threshold ≥ 0.95 indicates samples with clonal (single
> strain) infections, while samples with an F<sub>WS</sub> \< 0.95 are
> considered highly likely to come from mixed strain infections,
> indicating within-host diversity.

Source : [Amegashie et
al. (2020)](https://doi.org/10.1186/s12936-020-03510-3).

The F<sub>WS</sub> must be calculated by population. So we split up the
VCF by country.

Version: vcftools v0.1.16,
[vcfdo](https://github.com/IDEELResearch/vcfdo).

``` bash
VCFIN=malariae_brasi.filtsnps.core.ploidy2.vcf.gz

while read -r POP
do
  # Calculate by pop
vcftools --gzvcf  $VCFIN --keep $POP.samples \
--recode --stdout | gzip -c > $POP.vcf.gz #create the file by country
echo "Sample  Fws Standard_errors nb_sites" > ./fws/fws_$POP.txt
vcfdo wsaf -i $POP.vcf.gz | vcfdo fws >> ./fws/fws_$POP.txt #calculate fws
done < poplist.txt

# For sample with one sample per population, just calculate the heterozigosity
cat ./pop1/*.samples > pop1.samples
vcftools --gzvcf $VCFIN --keep pop1.samples \
--het --out ./fws/het_pop1.txt
```

We remove the individuals with a F<sub>WS</sub> or heterozigosity \<=
0.85 : 214 individuals remained.

## Remove related samples (IBD)

Version: bcftools v1.16, [hmmIBD](https://github.com/glipsnort/hmmIBD),
R v4.5.2, python v3.12.

> Highly related samples and clones can generate spurious signals of
> population structure, bias estimators of population genetic variation,
> and violate the assumptions of the model-based population genetic
> approaches ([Wang 2018](https://doi.org/10.1111/1755-0998.12708)). The
> relatedness between haploid genotype pairs was measured by estimating
> the pairwise fraction of the genome identical by descent (*IBD*)
> between strains within populations.

The IBD must be calculated by countries (only with n\>=2).

``` bash
# Transform in hmmIBD format
bcftools view malariae_brasi.filtsnps.core.ploidy2.vcf.gz \
-S keep_FWSHet.samples -o malariae_brasi.fws.ploidy2.vcf.gz

python vcf2hmm.py malariae_brasi.fws.ploidy2.vcf.gz malariae_brasi

# Running hmmIBD
hmmIBD -i malariae_brasi_seq.txt -o IBD_malariae_brasi -f malariae_brasi_freq.txt
```

The script `vcf2hmm.py` is from [hmmIBD
GitHub](https://github.com/glipsnort/hmmIBD).

Read the data and filter it

``` r
library(readr)
# Read the data from hmmIBD
IBD_pairs <-read_delim("./Data/IBD/IBD_malariae_brasi.hmm_fract.txt", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)

# Add the informations 
bam_paths <- read_delim("./Data/nofilt_meta.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
bam_paths$pop<-paste0(bam_paths$Species,"_",bam_paths$Country)
bam_paths$pop<-gsub(" ","",bam_paths$pop)
bam_paths$pop<-gsub("P\\.","",bam_paths$pop)

# keep only the pairs from the same populations
library(dplyr)
info_a<-bam_paths %>% rename_with(~paste0(., "_a"), grep("^[A-Z]*", names(.)))
IBD_pairs<-inner_join(IBD_pairs, info_a, join_by(sample1==`Other name_a`))

info_b<-bam_paths %>% rename_with(~paste0(., "_b"), grep("^[A-Z]*", names(.)))
IBD_pairs<-inner_join(IBD_pairs, info_b, join_by(sample2==`Other name_b`))

IBD_all<-subset(IBD_pairs, pop_a==pop_b)

### Filter
total_list<-unique(c(IBD_all$sample1,IBD_all$sample2)) # List of all the samples
# Only keep pair of individuals with IBD>0.5
fam_IBD<-subset(IBD_all, IBD_all$fract_sites_IBD>0.5)

#Assign family factor by individuals
clst = data.frame(ind = c(as.character(fam_IBD$sample1[1]), as.character(fam_IBD$sample2[1])), grp = c(1,1)) # initialize data.frame
for(i in 2:dim(fam_IBD)[1]){
  if(length(which(as.character(fam_IBD$sample1[i])==clst$ind))>0){
    tmp = data.frame(ind = c(as.character(fam_IBD$sample1[i]), as.character(fam_IBD$sample2[i])), grp = c(clst$grp[which(as.character(fam_IBD$sample1[i])==clst$ind)],clst$grp[which(as.character(fam_IBD$sample1[i])==clst$ind)]))
    clst = rbind(clst, tmp)
  } else if(length(which(as.character(fam_IBD$sample2[i])==clst$ind))>0){
    tmp = data.frame(ind = c(as.character(fam_IBD$sample1[i]), as.character(fam_IBD$sample2[i])), grp = c(clst$grp[which(as.character(fam_IBD$sample2[i])==clst$ind)],clst$grp[which(as.character(fam_IBD$sample2[i])==clst$ind)]))
    clst = rbind(clst, tmp)
  } else {
    tmp = data.frame(ind = c(as.character(fam_IBD$sample1[i]), as.character(fam_IBD$sample2[i])), grp = c(max(clst$grp)+1,max(clst$grp)+1))
    clst = rbind(clst, tmp)
  }
  clst = unique(clst)
}

# add missing data information
library(readr)
median_cov <- read_delim("./Data/malariae_brasi.info.ploidy2.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1, show_col_types = FALSE)

data_fam<-inner_join(clst, median_cov)

#keep the individual in each family with the higher median coverage
unrelated<-data_fam %>% 
    group_by(grp) %>% 
    slice(which.min(fmiss))

# make the dataframe for the samples discarded
related <- data_fam[!data_fam$ind %in% unrelated$ind,]
```

Filter the sample according to the median coverage. For each pair of
inbred samples (IBD \> 0.5), keep the sample with the higher median
coverage. 183 samples remain.
