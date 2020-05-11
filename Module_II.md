## Module II.  Imputation

#### Information
* Created by: GP2 Training and Networking

## Table of Contents
#### [0. Getting Started](#0)
#### [1. Post-QC data formatting](#1)
#### [1a. Check data with imputation reference SNP list](#2)
#### [1b. Convert to VCF](#3)
#### [1c. Sort and compress VCF](#4)
#### [1d. CheckVCF.py (optional)](#5)
#### [2. Generating Softcall + Hardcall](#6)

---

<a id="0"></a>
## 0. Getting Started

Welcome to module II of GP2 Bioinformatics Course. Please note that this markdown notebook is designed to be followed along with the video lesson.

Before you begin, please make sure that the following programs and scripts are installed or downloaded:

*Programs*
* Bash shell
* PLINK 1.9
* BCFTools
* 7zip or other decompression software
* Perl
* Python 2.7 (optional)
* R 3.6 (optional)

*Scripts*
* HRC-1000G-check-bim-v4.2.13.zip (unzipped) by William Rayner: https://www.well.ox.ac.uk/~wrayner/tools/
* HRC Panel (HRC.r1-1.GRCh37.wgs.mac5.sites.tab) http://www.haplotype-reference-consortium.org/site 
* CheckVCF (optional) https://github.com/zhanxw/checkVCF

As well as quality controlled genotype data in plink binary. For further information on quality controlled plink binary file, please check out Module I.

Our workspace folder structure is as follows:

    GP2_imputation_module
    ├── FILTERED_demofile
    │   ├── FILTERED.test.bed
    │   ├── FILTERED.test.bim
    │   ├── FILTERED.test.fam
    │   └── FILTERED.test.log
    └── Workspace
        ├── checkVCF.py
        ├── filter_variants_maf_softcall_r2.R
        ├── HRC-1000G-check-bim.pl
        ├── HRC.r1-1.GRCh37.wgs.mac5.sites.tab
        └── GP2_module2.ipynb

Your folder or workspace does not have to look like this. It's just given to you so you can follow along the demo.

Our starting input data are PLINK binary files `FILTERED_demofile/FILTERED.test{.bed,bim,fam}`.

<a id="1"></a>
## 1. Post-QC data formatting

This should start with data that's already been quality controlled as described by module 1, including but not limited to population structure and Hardy-Weinberg.

**Four steps**:

* Check data with imputation reference SNP list
* Convert to VCF
* Compress and sort VCF
* CheckVCF.py for final QC (optional)

<a id="2"></a>
### 1a. Check data with imputation reference SNP list

Using William Rayner's tool: https://www.well.ox.ac.uk/~wrayner/tools/

Reference list we will use it against: http://www.haplotype-reference-consortium.org/site

(see the .tab file)


```bash
FILENAME=FILTERED.test
```


```bash
pwd
```

    /Users/kimjoj/Documents/GP2_imputation_module/workspace



```bash
export PATH=$PATH:/Users/kimjoj/bash_programs/plink:/Users/kimjoj/bash_programs/bcftools
```


```bash
echo $PATH
```

    /Users/kimjoj/anaconda3/bin:/Users/kimjoj/anaconda3/bin:/Users/kimjoj/anaconda3/condabin:/usr/bin:/bin:/usr/sbin:/sbin:/Users/kimjoj/bash_programs/plink:/Users/kimjoj/bash_programs/bcftools



```bash
# William Rayner script requires frequency output as an input
plink --bfile ../FILTERED_demofile/$FILENAME --freq --out $FILENAME

perl HRC-1000G-check-bim.pl -b ../FILTERED_demofile/$FILENAME.bim -f $FILENAME.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
```

    PLINK v1.90b6.17 64-bit (28 Apr 2020)          www.cog-genomics.org/plink/1.9/
    (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to FILTERED.test.log.
    Options in effect:
      --bfile ../FILTERED_demofile/FILTERED.test
      --freq
      --out FILTERED.test
    
    16384 MB RAM detected; reserving 8192 MB for main workspace.
    839 variants loaded from .bim file.
    2000 people (1040 males, 960 females) loaded from .fam.
    2000 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 2000 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    --freq: Allele frequencies (founders only) written to FILTERED.test.frq .
    
    
             Script to check plink .bim files against HRC/1000G for
            strand, id names, positions, alleles, ref/alt assignment
                        William Rayner 2015-2020
                            wrayner@well.ox.ac.uk
    
                                Version 4.2.13
    
    
    Options Set:
    Reference Panel:             HRC
    Bim filename:                ../FILTERED_demofile/FILTERED.test.bim
    Reference filename:          HRC.r1-1.GRCh37.wgs.mac5.sites.tab
    Allele frequencies filename: FILTERED.test.frq
    Plink executable to use:     plink
    
    Chromosome flag set:         No
    Allele frequency threshold:  0.2
    
    Path to plink bim file: /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile
    
    Reading HRC.r1-1.GRCh37.wgs.mac5.sites.tab
    2000000
    4000000
    6000000
    8000000
    10000000
    12000000
    14000000
    16000000
    18000000
    20000000
    22000000
    24000000
    26000000
    28000000
    30000000
    32000000
    34000000
    36000000
    38000000
    40000000
    40405506
     ...Done
    
    
    Details written to log file: /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/LOG-FILTERED.test-HRC.txt
    
    Creating variant lists
    /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/Force-Allele1-FILTERED.test-HRC.txt
    /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/Strand-Flip-FILTERED.test-HRC.txt
    /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/ID-FILTERED.test-HRC.txt
    /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/Position-FILTERED.test-HRC.txt
    /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/Chromosome-FILTERED.test-HRC.txt
    /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/Exclude-FILTERED.test-HRC.txt
    /Users/kimjoj/Documents/GP2_imputation_module/workspace/../FILTERED_demofile/FreqPlot-FILTERED.test-HRC.txt
    
    
    Matching to HRC
    
    Position Matches
     ID matches HRC 0
     ID Doesn't match HRC 839
     Total Position Matches 839
    ID Match
     Position different from HRC 0
    No Match to HRC 0
    Skipped (MT) 0
    Total in bim file 839
    Total processed 839
    
    Indels 0
    
    SNPs not changed 0
    SNPs to change ref alt 839
    Strand ok 839
    Total Strand ok 839
    
    Strand to change 0
    Total checked 839
    Total checked Strand 839
    Total removed for allele Frequency diff > 0.2 0
    Palindromic SNPs with Freq > 0.4 0
    
    
    Non Matching alleles 0
    ID and allele mismatching 0; where HRC is . 0
    Duplicates removed 0
    
    
    Writing plink commands to: Run-plink.sh



```bash
# Creates Run-plink.sh
# Check log file to check for potential issues

sh ../FILTERED_demofile/Run-plink.sh
```

Now we should see `FILTERED.test-updated-chr` for each chromosomes.

<a id="3"></a>
### 1b. Convert to VCF

When converting to VCF, we will split this to different chromosomes as required by Michigan Imputation Server (it also makes files smaller and easier to handle)


```bash
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  plink --bfile ../FILTERED_demofile/$FILENAME-updated-chr$chnum --recode vcf --chr $chnum --out $FILENAME$chnum
done
```

<a id="4"></a>
### 1c. Sort and compress VCF

MIS also requires bgzip compression.


```bash
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  bcftools sort $FILENAME$chnum.vcf -Oz -o pre_impute_$FILENAME$chnum.vcf.gz
  #vcf-sort $FILENAME$chnum.vcf | bgzip -c >  pre_impute_$FILENAME$chnum.vcf.gz
done
```

    Writing to /tmp/bcftools-sort.h7vLd2
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.XdDazc
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.CeHalg
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.lNQZYE
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.ryX7b2
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.txlG2o
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.yhxTae
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.V0RYr2
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.m0A9FM
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.aciYF4
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.U4KQyQ
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.GpfFKE
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.DzPRS5
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.hjGSHW
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.XvMhft
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.Qxq69l
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.K0S5XZ
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.rGgVSo
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.quusOg
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.1LYOb8
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.9QKScU
    Merging 1 temporary files
    Cleaning
    Done
    Writing to /tmp/bcftools-sort.rWZNXw
    Merging 1 temporary files
    Cleaning
    Done


 <a id="5"></a>
### 1d. CheckVCF.py (optional)

Found here: https://github.com/zhanxw/checkVCF

Additional validation checks such as indel sites and ref allele consistency


```bash
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  python2.7 checkVCF.py -r hs37d5.fa -o $FILENAME$chnum pre_impute_$FILENAME$chnum.vcf.gz
done
```

    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test1.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 61 ] lines processed
    Examine [ 8 ] VCF header lines, [ 53 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test1.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test1.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test1.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test1.check.af ]
    [ 13 ] Monomorphic sites are outputted to [ FILTERED.test1.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test1.check.log FILTERED.test1.check.dup FILTERED.test1.check.noSnp FILTERED.test1.check.ref FILTERED.test1.check.geno FILTERED.test1.check.af FILTERED.test1.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test2.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 73 ] lines processed
    Examine [ 8 ] VCF header lines, [ 65 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test2.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test2.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test2.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test2.check.af ]
    [ 16 ] Monomorphic sites are outputted to [ FILTERED.test2.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test2.check.log FILTERED.test2.check.dup FILTERED.test2.check.noSnp FILTERED.test2.check.ref FILTERED.test2.check.geno FILTERED.test2.check.af FILTERED.test2.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test3.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 70 ] lines processed
    Examine [ 8 ] VCF header lines, [ 62 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test3.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test3.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test3.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test3.check.af ]
    [ 19 ] Monomorphic sites are outputted to [ FILTERED.test3.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test3.check.log FILTERED.test3.check.dup FILTERED.test3.check.noSnp FILTERED.test3.check.ref FILTERED.test3.check.geno FILTERED.test3.check.af FILTERED.test3.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test4.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 40 ] lines processed
    Examine [ 8 ] VCF header lines, [ 32 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test4.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test4.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test4.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test4.check.af ]
    [ 11 ] Monomorphic sites are outputted to [ FILTERED.test4.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test4.check.log FILTERED.test4.check.dup FILTERED.test4.check.noSnp FILTERED.test4.check.ref FILTERED.test4.check.geno FILTERED.test4.check.af FILTERED.test4.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test5.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 102 ] lines processed
    Examine [ 8 ] VCF header lines, [ 94 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test5.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test5.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test5.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test5.check.af ]
    [ 31 ] Monomorphic sites are outputted to [ FILTERED.test5.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test5.check.log FILTERED.test5.check.dup FILTERED.test5.check.noSnp FILTERED.test5.check.ref FILTERED.test5.check.geno FILTERED.test5.check.af FILTERED.test5.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test6.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 83 ] lines processed
    Examine [ 8 ] VCF header lines, [ 75 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test6.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test6.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test6.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test6.check.af ]
    [ 26 ] Monomorphic sites are outputted to [ FILTERED.test6.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test6.check.log FILTERED.test6.check.dup FILTERED.test6.check.noSnp FILTERED.test6.check.ref FILTERED.test6.check.geno FILTERED.test6.check.af FILTERED.test6.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test7.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 54 ] lines processed
    Examine [ 8 ] VCF header lines, [ 46 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test7.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test7.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test7.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test7.check.af ]
    [ 11 ] Monomorphic sites are outputted to [ FILTERED.test7.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test7.check.log FILTERED.test7.check.dup FILTERED.test7.check.noSnp FILTERED.test7.check.ref FILTERED.test7.check.geno FILTERED.test7.check.af FILTERED.test7.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test8.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 19 ] lines processed
    Examine [ 8 ] VCF header lines, [ 11 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test8.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test8.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test8.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test8.check.af ]
    [ 2 ] Monomorphic sites are outputted to [ FILTERED.test8.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test8.check.log FILTERED.test8.check.dup FILTERED.test8.check.noSnp FILTERED.test8.check.ref FILTERED.test8.check.geno FILTERED.test8.check.af FILTERED.test8.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test9.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 51 ] lines processed
    Examine [ 8 ] VCF header lines, [ 43 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test9.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test9.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test9.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test9.check.af ]
    [ 16 ] Monomorphic sites are outputted to [ FILTERED.test9.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test9.check.log FILTERED.test9.check.dup FILTERED.test9.check.noSnp FILTERED.test9.check.ref FILTERED.test9.check.geno FILTERED.test9.check.af FILTERED.test9.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test10.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 26 ] lines processed
    Examine [ 8 ] VCF header lines, [ 18 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test10.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test10.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test10.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test10.check.af ]
    [ 4 ] Monomorphic sites are outputted to [ FILTERED.test10.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test10.check.log FILTERED.test10.check.dup FILTERED.test10.check.noSnp FILTERED.test10.check.ref FILTERED.test10.check.geno FILTERED.test10.check.af FILTERED.test10.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test11.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 25 ] lines processed
    Examine [ 8 ] VCF header lines, [ 17 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test11.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test11.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test11.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test11.check.af ]
    [ 3 ] Monomorphic sites are outputted to [ FILTERED.test11.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test11.check.log FILTERED.test11.check.dup FILTERED.test11.check.noSnp FILTERED.test11.check.ref FILTERED.test11.check.geno FILTERED.test11.check.af FILTERED.test11.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test12.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 41 ] lines processed
    Examine [ 8 ] VCF header lines, [ 33 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test12.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test12.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test12.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test12.check.af ]
    [ 14 ] Monomorphic sites are outputted to [ FILTERED.test12.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test12.check.log FILTERED.test12.check.dup FILTERED.test12.check.noSnp FILTERED.test12.check.ref FILTERED.test12.check.geno FILTERED.test12.check.af FILTERED.test12.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test13.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 30 ] lines processed
    Examine [ 8 ] VCF header lines, [ 22 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test13.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test13.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test13.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test13.check.af ]
    [ 3 ] Monomorphic sites are outputted to [ FILTERED.test13.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test13.check.log FILTERED.test13.check.dup FILTERED.test13.check.noSnp FILTERED.test13.check.ref FILTERED.test13.check.geno FILTERED.test13.check.af FILTERED.test13.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test14.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 29 ] lines processed
    Examine [ 8 ] VCF header lines, [ 21 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test14.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test14.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test14.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test14.check.af ]
    [ 4 ] Monomorphic sites are outputted to [ FILTERED.test14.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test14.check.log FILTERED.test14.check.dup FILTERED.test14.check.noSnp FILTERED.test14.check.ref FILTERED.test14.check.geno FILTERED.test14.check.af FILTERED.test14.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test15.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 41 ] lines processed
    Examine [ 8 ] VCF header lines, [ 33 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test15.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test15.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test15.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test15.check.af ]
    [ 12 ] Monomorphic sites are outputted to [ FILTERED.test15.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test15.check.log FILTERED.test15.check.dup FILTERED.test15.check.noSnp FILTERED.test15.check.ref FILTERED.test15.check.geno FILTERED.test15.check.af FILTERED.test15.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test16.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 57 ] lines processed
    Examine [ 8 ] VCF header lines, [ 49 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test16.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test16.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test16.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test16.check.af ]
    [ 13 ] Monomorphic sites are outputted to [ FILTERED.test16.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test16.check.log FILTERED.test16.check.dup FILTERED.test16.check.noSnp FILTERED.test16.check.ref FILTERED.test16.check.geno FILTERED.test16.check.af FILTERED.test16.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test17.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 50 ] lines processed
    Examine [ 8 ] VCF header lines, [ 42 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test17.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test17.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test17.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test17.check.af ]
    [ 6 ] Monomorphic sites are outputted to [ FILTERED.test17.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test17.check.log FILTERED.test17.check.dup FILTERED.test17.check.noSnp FILTERED.test17.check.ref FILTERED.test17.check.geno FILTERED.test17.check.af FILTERED.test17.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test18.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 10 ] lines processed
    Examine [ 8 ] VCF header lines, [ 2 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test18.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test18.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test18.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test18.check.af ]
    [ 0 ] Monomorphic sites are outputted to [ FILTERED.test18.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test18.check.log FILTERED.test18.check.dup FILTERED.test18.check.noSnp FILTERED.test18.check.ref FILTERED.test18.check.geno FILTERED.test18.check.af FILTERED.test18.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test19.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 31 ] lines processed
    Examine [ 8 ] VCF header lines, [ 23 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test19.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test19.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test19.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test19.check.af ]
    [ 7 ] Monomorphic sites are outputted to [ FILTERED.test19.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test19.check.log FILTERED.test19.check.dup FILTERED.test19.check.noSnp FILTERED.test19.check.ref FILTERED.test19.check.geno FILTERED.test19.check.af FILTERED.test19.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test20.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 80 ] lines processed
    Examine [ 8 ] VCF header lines, [ 72 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test20.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test20.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test20.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test20.check.af ]
    [ 19 ] Monomorphic sites are outputted to [ FILTERED.test20.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test20.check.log FILTERED.test20.check.dup FILTERED.test20.check.noSnp FILTERED.test20.check.ref FILTERED.test20.check.geno FILTERED.test20.check.af FILTERED.test20.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test21.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 19 ] lines processed
    Examine [ 8 ] VCF header lines, [ 11 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test21.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test21.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test21.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test21.check.af ]
    [ 4 ] Monomorphic sites are outputted to [ FILTERED.test21.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test21.check.log FILTERED.test21.check.dup FILTERED.test21.check.noSnp FILTERED.test21.check.ref FILTERED.test21.check.geno FILTERED.test21.check.af FILTERED.test21.check.mono
    checkVCF.py -- check validity of VCF file for meta-analysis
    version 1.4 (20140115)
    contact zhanxw@umich.edu or dajiang@umich.edu for problems.
    Python version is [ 2.7.16.final.0 ] 
    Begin checking vcfFile [ pre_impute_FILTERED.test22.vcf.gz ]
    ---------------     REPORT     ---------------
    Total [ 23 ] lines processed
    Examine [ 8 ] VCF header lines, [ 15 ] variant sites, [ 2000 ] samples
    [ 0 ] duplicated sites
    [ 0 ] NonSNP site are outputted to [ FILTERED.test22.check.nonSnp ]
    [ 0 ] Inconsistent reference sites are outputted to [ FILTERED.test22.check.ref ]
    [ 0 ] Variant sites with invalid genotypes are outputted to [ FILTERED.test22.check.geno ]
    [ 0 ] Alternative allele frequency > 0.5 sites are outputted to [ FILTERED.test22.check.af ]
    [ 6 ] Monomorphic sites are outputted to [ FILTERED.test22.check.mono ]
    ---------------     ACTION ITEM     ---------------
    * No error found by checkVCF.py, thank you for cleanning VCF file.
    * Upload these files to the ftp server (so we can double check): FILTERED.test22.check.log FILTERED.test22.check.dup FILTERED.test22.check.noSnp FILTERED.test22.check.ref FILTERED.test22.check.geno FILTERED.test22.check.af FILTERED.test22.check.mono


Looks like everything worked out. Let's move the results to somewhere easier to find.


```bash
mv *.vcf.gz ../postformat/
```

 <a id="6"></a>
## 2. Generating Softcall + hardcall binaries

After Michigan Imputation Server, we have our imputed data, but not imputed sites are equal... We need to determine (arbitrary) cutoff points for both minor allele frequency and Rsq.

Softcall: rsq > 0.3

Hardcall: rsq > 0.8

---
### Optional Step: what sites pass our MAF and Rsq cutoffs?

Copy and save the following as `filter_variants_maf_r2_updated.R`

    if (!"tidyverse" %in% rownames(installed.packages())) {
      install.packages("tidyverse", repos = "https://cloud.r-project.org")
    }
    if (!"data.table" %in% rownames(installed.packages())) {
      install.packages("data.table", repos = "https://cloud.r-project.org")
    }

    library(tidyverse)
    library(data.table)

    arg <- commandArgs(trailingOnly = T)

    # define function we will use
    filter_info <- function (chr, maf, rsq) {
      # Read the info file of specific chromosome
      data <- paste0("chr", chr,".info") %>% fread()
      # Filter it to maf and rsq
      dat <- data[data$MAF >= maf & data$Rsq >= rsq]
      # Generate chromosome, bp, and range from the "SNP" column
      dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
      dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
      dat$range <- paste0(dat$chr, ":", dat$bp, "-", dat$bp)
      # writing files
      dat[,c("SNP","ALT_Frq","Rsq")] %>% fwrite(
        paste0("maf", gsub("0\\.", "", maf), "rsq", gsub("0\\.", "", rsq), "minimums_chr", chr, ".info"),
        row.names = F, quote = F, sep = "\t"
      )
      dat[,c("range")] %>% fwrite(
        paste0("maf", gsub("0\\.", "", maf), "rsq", gsub("0\\.", "", rsq), "minimums_chr", chr, ".txt"),
        col.names = F, row.names = F, quote = F, sep = "\t"
      )
      # report number of SNPs that pass the filter
      paste0("Number of SNPs in chromosome ", chr, " after filter (MAF >= " , maf, ", Rsq >= ", rsq, "): ", nrow(dat)) %>% print()
    }

    # SET MAF and RSQ filters here 

    lapply(1:22, filter_info, maf = arg[1], rsq = arg[2])


```bash
gunzip *.info.gz

# SYNTAX: Rscript filter_variants_maf_r2_updated.R [MAF cutoff] [Rsq cutoff]
# following will use MAF cutoff of 0.001 and Rsq cutoff of 0.8
Rscript filter_variants_maf_r2_updated.R 0.001 0.8
```

Produces two outputs for each chromosomes. They contain SNPs that pass the given MAF and Rsq filters.
1. .info file contains SNP name (chr:bp:ref:alt), MAF, and Rsq
2. .txt file contains the chr:bp range of the individual SNPs


```bash
cat maf001rsq03minimums_chr1.info
```

    SNP	ALT_Frq	Rsq
    1:21012549:G:A	0.00195	0.97271
    1:21012595:C:T	0.00117	0.93975
    1:21012609:G:A	0.00519	0.98883
    1:27688743:T:C	0.00472	0.98031
    1:154909820:A:G	0.0012	0.95889
    1:154966651:G:A	0.00166	0.9461
    1:169953815:A:C	0.00271	0.98424
    1:175372714:T:G	0.00346	0.98766
    1:175375355:T:C	0.00349	0.97695
    1:205595114:T:C	0.00221	0.98292
    1:205714901:T:C	0.00197	0.98727



```bash
cat maf001rsq03minimums_chr1.txt
```

    1:21012549-21012549
    1:21012595-21012595
    1:21012609-21012609
    1:27688743-27688743
    1:154909820-154909820
    1:154966651-154966651
    1:169953815-169953815
    1:175372714-175372714
    1:175375355-175375355
    1:205595114-205595114
    1:205714901-205714901


---
Now let's generate a softcall. The script below converts the `dose.vcf.gz` files to Plink binaries while filtering by given Rsq cutoff value.


```bash
# make the output folder
mkdir plink_files_soft

#change "--qual-threshold" to desired rsq threshold. in this case, it is 0.3 because we are generating softcall

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  plink --vcf chr$chnum.dose.vcf.gz --qual-scores chr$chnum.info 7 1 1 --qual-threshold 0.3 --make-bed --out plink_files_soft_2/chr$chnum  --double-id
done
```

Change the number that comes after `--qual-threshold` to change the Rsq cutoff value. In this case since it is softcalls, it is set to `0.3`. If you want it to be more stringent, say a hardcall, then set it to `0.8`.


```bash
cd plink_files_soft
ls | grep ".bim" > merge_list.txt
sed -i 's/.bim//g' merge_list.txt
plink --merge-list merge_list.txt --make-bed --out IMPUTED.SOFTCALLS.Demo
rm chr*.bim
rm chr*.log
rm chr*.fam
rm chr*.bed
rm chr*.nosex
cd ..
```

This merges all Plink binaries separated by chrosomes into one single binary dataset.


```bash
# recompress .info file
gzip *.info
```
