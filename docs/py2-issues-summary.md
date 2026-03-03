# LDSC Python 2 Repo Issues — Adversarial Summary (builk/ldsc)

This document summarizes the most discussed issues in the original Python LDSC repository (buliK/ldsc) using a simple heuristic based primarily on **comment count**. Pull requests are excluded.

Data source: GitHub issues list (state=all). Total issues scanned: 449.

## Heuristic
- Primary: number of comments (descending).
- Secondary: reactions total (when present).

## Top Issues by Discussion
### 1. #26 — munge_sumstats.py
- State: closed
- Comments: 45
- Created: 2015-02-20T13:49:50Z
- Updated: 2016-02-19T19:24:44Z
- Excerpt: Hello, Thank you for making this interesting piece of software I'm keen on using it on my own data sets. I've ran into a spot of bother using the munge_sumstats.py provided (downloaded on the 20/01/2015). I'm trying t...

### 2. #369 — Getting module parse has no attibute error
- State: open
- Comments: 31
- Created: 2023-02-24T18:45:37Z
- Updated: 2024-12-06T01:57:28Z
- Excerpt: Hello, I am trying to calcualate heritbaility and I am getting this issue. May I know how to resolve this? Code: python3 dsc.py --h2 sbp_hapmap.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out sbp...

### 3. #374 — Why do I keep reporting an error: Could not find SNP column？
- State: open
- Comments: 25
- Created: 2023-03-24T03:33:10Z
- Updated: 2024-12-03T10:38:28Z
- Excerpt: (ldsc) python munge_sumstats.py \ > --sumstats /mnt/ndisk1/Student/zhouyi/GWAS_tools/TWAS/data/wellbeing_sum/WBresults.txt \ > --N 491455 \ > --signed-sumstats beta,0 \ > --snp snpid \ > --out SWB \ > --merge-alleles ...

### 4. #43 — fatal error while running munge_sumstats.py
- State: open
- Comments: 24
- Created: 2016-02-26T19:17:03Z
- Updated: 2025-02-04T03:46:43Z
- Excerpt: Hello, After successful installation of ldsc on my system, I tried to run munge_sumstats.py script to format my input files to ldsc format files. I thought this will be a smooth process but it's not. While running mun...

### 5. #66 — Error converting summary statistics
- State: open
- Comments: 20
- Created: 2016-12-15T18:08:28Z
- Updated: 2023-11-07T07:25:41Z
- Excerpt: Getting the following error. Any help would be appreciated! Thanks. ------------------------------- Call: ./munge_sumstats.py \ --out outfile\ --merge-alleles w_hm3.snplist \ --N 4810.0 \ --sumstats infile.txt \ --ign...

### 6. #425 — Is it necessary to use the --merge-alleles parameter when doing genetic correlation in ldsc?
- State: open
- Comments: 16
- Created: 2024-03-21T05:26:52Z
- Updated: 2025-12-09T19:50:32Z
- Excerpt: Hi, I have two sumstats, before doing genetic correlation, I have intersected the SNPs of the two data and generated the sumstats.gz files. When I do rg, it reports this error: ./ldsc.py \ --ref-ld-chr eur_w_ld_chr/ \...

### 7. #438 — Prop._h2 is negative
- State: open
- Comments: 15
- Created: 2024-06-21T10:45:29Z
- Updated: 2024-07-29T16:12:25Z
- Excerpt: Hi, When I did the prtitioned heritability with 1000G_Phase3_baselineLD_v2.2_ldscores.tgz, the Prop._h2 was negative and it was very significant. Is this a bug?

### 8. #37 — median value issue for munge_sumstats.py
- State: closed
- Comments: 14
- Created: 2015-10-22T14:42:30Z
- Updated: 2026-02-26T19:53:51Z
- Excerpt: I have problem using munge_sumstats.py for one of gwas data, I used the effect size as the summary statistics, the errors shows below: Traceback (most recent call last): File "/nas02/home/k/x/kxia/software/ldsc/munge_...

### 9. #145 — FIXED - Fail to converst summary statistic in .sumstats format: munge_sumstats is taking hours 
- State: open
- Comments: 13
- Created: 2019-02-21T08:26:40Z
- Updated: 2025-08-19T07:41:35Z
- Excerpt: I am trying to reproduce the example provided in: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation In particular, I downloaded both the summary statistics file: wget www.med.unc.edu/pgc/files/re...

### 10. #366 — IndexError: list index out of range
- State: open
- Comments: 13
- Created: 2022-12-05T03:34:03Z
- Updated: 2025-11-03T20:12:52Z
- Excerpt: Getting the following error. Any help would be appreciated! Thanks. * LD Score Regression (LDSC) * Version 1.0.1 * (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane * Broad Institute of MIT and Harvard / MIT De...

### 11. #138 — Prop._h2 always is 1.0 use my own ldscores based on my own annot
- State: open
- Comments: 12
- Created: 2019-01-02T18:38:04Z
- Updated: 2021-06-17T13:58:05Z
- Excerpt: How to explain my result: Prop._h2 = 1.0 using my ldscore based on my own annot? Many thanks.

### 12. #371 — link to reference files does not work
- State: open
- Comments: 11
- Created: 2023-03-10T04:39:19Z
- Updated: 2024-06-04T14:34:20Z
- Excerpt: The link to w_hm3.snplist as well as EUR ld reference files do not seem to work https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist....

### 13. #39 — ./ldsc.py -h ERROR
- State: closed
- Comments: 11
- Created: 2015-10-23T15:34:58Z
- Updated: 2015-11-20T09:36:45Z
- Excerpt: Hi, I just installed ldsc, seem to have all the python packages installed, but still get this message: -bash-4.1$ ./ldsc.py -h Traceback (most recent call last): File "./ldsc.py", line 13, in <module> import ldscore.p...

### 14. #158 — LinAlgError: Singular matrix
- State: open
- Comments: 10
- Created: 2019-07-08T16:03:24Z
- Updated: 2022-12-19T13:35:27Z
- Excerpt: Hello, I tried to use the LDSC for Partitioned Heritability from Continuous Annotations, but we got the following error information. Can you help to check what's the reason for this kind of error? Thanks. ************...

### 15. #327 — munge_sumstats.py_error "ValueError: No object to concatenate"
- State: open
- Comments: 10
- Created: 2021-11-03T19:18:06Z
- Updated: 2025-03-02T06:36:03Z
- Excerpt: I am trying to to cross-trait LDSC for a multivariate analysis. In trying to use the munge_sumstats.py to format the statistics, I am getting this error on multiple files. Any suggestions? ****************************...

### 16. #407 — Suitable w_hm3.snplist for AFR?
- State: open
- Comments: 10
- Created: 2023-10-11T17:17:20Z
- Updated: 2025-07-08T16:47:59Z
- Excerpt: Hi! I've come across a LD reference panel for AFR. But I need an appropriate SNP list for --merge-alleles. Does anyone know how to generate something suitable? Thank you!

### 17. #309 — error
- State: open
- Comments: 10
- Created: 2021-09-10T07:11:53Z
- Updated: 2025-04-14T08:26:12Z
- Excerpt: Call: ./ldsc.py \ --ref-ld-chr /home/huawei/xb/software/ldsc/eur_w_ld_chr \ --out test \ --rg A1.sumstats.gz,phylum.B1.sumstats.gz \ --w-ld-chr /home/huawei/xb/software/ldsc/eur_w_ld_chr Beginning analysis at Fri Sep ...

### 18. #155 — negative h2 and Ratio > 1
- State: open
- Comments: 10
- Created: 2019-05-08T04:56:33Z
- Updated: 2022-06-27T22:59:41Z
- Excerpt: Dear Brendan, I just ran your code on my example and had the wired result as below: Total Observed scale h2: -0.1825 (0.0625) Lambda GC: 2.5879 Mean Chi^2: 3.5595 Intercept: 3.8665 (0.052) Ratio: 1.1199 (0.0203) The h...

### 19. #147 — Cell type analysis: IOError: Could not open Cahoy.1.l2.ldscore[./gz/bz2]
- State: open
- Comments: 10
- Created: 2019-03-06T20:11:02Z
- Updated: 2025-08-23T04:40:45Z
- Excerpt: I am trying to do the cell type analyses following the example with the Cahoy cell types and consistently get the error that Cahoy.1.l2.ldscore[./gz/bz2] cannot be opened. All files are executable, I tried moving the ...

### 20. #117 — LD score regression trait correlation ERROR:invalid value encountered in sqrt
- State: open
- Comments: 10
- Created: 2018-05-14T16:32:29Z
- Updated: 2018-05-18T21:26:16Z
- Excerpt: Hi, I am running LDSR correlation between two traits (t1d and w6), and when I run the analysis considering an intercept I get the following error:FloatingPointError: invalid value encountered in sqrt Since there is no...

## Label Frequency in Top Issues
- (No labels on top issues)

## Notes
- This is not an exhaustive thematic analysis; it is a triage view to show what users discuss most.
- If you want a deeper thematic grouping (e.g., installation, M issues, allele matching, rg stability), I can expand this.
