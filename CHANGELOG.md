# Changelog

## [0.6.1](https://github.com/RIVM-bioinformatics/SARS2seq/compare/v0.6.0...v0.6.1) (2022-08-15)


### Bug Fixes

* write AmpliGone output to logfile instead of stdout ([ff16fe4](https://github.com/RIVM-bioinformatics/SARS2seq/commit/ff16fe4f6436e56d31410bac45ad13b4d6fefef3))


### Dependencies

* simplify the `QC_and_cleanup` environment ([8a8059e](https://github.com/RIVM-bioinformatics/SARS2seq/commit/8a8059ee2131988bb6432bb5fa4dbe5d823d239e))
* update Nextclade to newer (lenient) version 2.4.x ([d5b90e2](https://github.com/RIVM-bioinformatics/SARS2seq/commit/d5b90e2d452f3999417fb655746f6768326e432a))
* Update required snakemake version to 7.x (pinned to 7.12.1) ([b171fc4](https://github.com/RIVM-bioinformatics/SARS2seq/commit/b171fc4c9dcc4882bb4e827ac73fd4b5c682389f))
* upgrade AmpliGone to version 1.1.0 ([ec005aa](https://github.com/RIVM-bioinformatics/SARS2seq/commit/ec005aa58f3d1a98ec45ed7f4ab669f27ef1a055))
* use the bioconda channel to install AmpliGone instead of pip ([8a8059e](https://github.com/RIVM-bioinformatics/SARS2seq/commit/8a8059ee2131988bb6432bb5fa4dbe5d823d239e))

## [0.6.0](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.5.4...v0.6.0) (2022-04-26)


### Features

* Update TrueConsense to 0.5.0 ([972c08a](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/972c08a3fd88e8ade9c6095ec84b3645306dcf05))


### Bug Fixes

* Update Ampligone to 1.0.1 ([cce0e31](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/cce0e31fa8b5a82a368936c526d2935d58620af6))

### [0.5.4](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.5.3...v0.5.4) (2022-04-08)


### Dependencies

* pin pangolin version number to 4.0.2 without auto-updating ([27e30d3](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/27e30d34940d852993f616b04feae231751bf607))

### [0.5.3](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.5.2...v0.5.3) (2022-04-05)


### Bug Fixes

* allow pangolin v4 empty output when there's low coverage in the typing sequence ([3237202](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/32372022b5f35a1d4332d4ad72eabc258b5ec2ea))

### [0.5.2](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.5.1...v0.5.2) (2022-04-04)


### Dependencies

* add python3.7 to Alignment env to stabilize environment during installation ([6f1387d](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/6f1387d9a10dbbcffe9dc664b7c160187288136b))
* change Typing rules and packages for compatibility with pangolin v4.* ([81a61d5](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/81a61d55ea320c516aaf1d6972638b22af8b3321))

### [0.5.1](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.5.0...v0.5.1) (2022-02-04)


### Bug Fixes

* finish an analysis when empty input fastq files are given ([7a30b18](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/7a30b18e1b149ac498cdd53902c9d1f2e258b1f7))


### Dependencies

* update AmpliGone 0.4.0 --> 0.4.1 ([cf64dc4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/cf64dc459be28e0b74abd1c2207d00e950b788c9))
* update AmpliGone 0.4.1 --> 0.4.3 ([28baf2a](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/28baf2a8becdc0e3945347a2483561fd119d728b))

## [0.5.0](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.4.2...v0.5.0) (2022-01-27)


### Features

* Write a short report with used/configured settings of the analysis ([ffe3141](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/ffe3141c2b3af151a108a2d768584e6d30cc1557))


### Bug Fixes

* properly exit with right exit codes ([4ec6307](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/4ec630723960f41f47c53e668d89f0f8a2cf1f8f))
* sucessfully finish an analysis when no primers are given ([4916731](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/491673192ea83293fb2fe639e833e7d712ea311e))


### Performance Improvements

* tweak read qc filtering settings ([27916d4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/27916d4a9dc6be9b402223eed36654886b55b541))


### Dependencies

* add FPDF version 1.7.2 to dependency list(s) ([1a192ec](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/1a192ecc146bd2d6eb1f07c0ea3f12331a042cf1))
* add urllib3 to dependency list for the auto-updater ([f3b9b53](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/f3b9b536228dd077a77e56d414b394c161ece147))
* include conda (lenient) version 4.11.x in main environment to circumvent python3.10 bug ([a7821f4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/a7821f463dd6157bea2e109f17cab9002cd138e1))
* remove pysamstats from 'clean' environment ([485c94e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/485c94ed732057a572f6a3e631aa5618f286e9cf))
* update (pangolin) snakemake 6.4.1 --> 6.8.0 ([1ffc79d](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/1ffc79dbb1669cf5df06b997ba28786690e9bde3))
* update AmpliGone version v0.3.3 --> v0.4.0 ([999b757](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/999b757518e40d6ea10fbd31ef35d7da81ff8772))
* update bcftools 1.12 --> 1.14 ([9108076](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/91080763fdfb65d4bb90372c7cef7d8302606f7b))
* update bedtools 2.29.2 --> 2.30.0 ([9108076](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/91080763fdfb65d4bb90372c7cef7d8302606f7b))
* update biopython 1.78 --> 1.79 ([9108076](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/91080763fdfb65d4bb90372c7cef7d8302606f7b))
* update biopython to version 1.79 in main environment ([43251f0](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/43251f09945d78cd3a7952d24256e8c703a915f9))
* update FastP to version 0.23.2 ([b134b64](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/b134b64644dbdb3af93c1bba3d23d3d50660e557))
* update fastqc 0.11.8 --> 0.11.9 ([485c94e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/485c94ed732057a572f6a3e631aa5618f286e9cf))
* update minimap2 & mappy 2.17 --> 2.24 ([cb99be4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/cb99be4361fbe232be214f0c1c556e560bf76201))
* update multiqc 1.9 --> 1.11 ([485c94e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/485c94ed732057a572f6a3e631aa5618f286e9cf))
* update nextclade 1.6.0 --> 1.9 + change to lenient versioning ([1ffc79d](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/1ffc79dbb1669cf5df06b997ba28786690e9bde3))
* update nextclade lenient version v1.9.x --> v1.10.x ([2e1cf9e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/2e1cf9ef2ef8c81bbb6219852df929a2f0fd684a))
* update pandas 1.2.3 --> 1.3.5 ([485c94e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/485c94ed732057a572f6a3e631aa5618f286e9cf))
* update pangolin 3.0 --> 3.1 + change to lenient versioning ([1ffc79d](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/1ffc79dbb1669cf5df06b997ba28786690e9bde3))
* update parmap 1.5.2 --> 1.5.3 ([485c94e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/485c94ed732057a572f6a3e631aa5618f286e9cf))
* update PyYaml 5.4.1 --> 6.0 and fix compatibility ([e9a1ef9](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/e9a1ef955948ec30ecb6960a69e9da7d4ef8bd27))
* update samtools 1.10 --> 1.14 ([cb99be4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/cb99be4361fbe232be214f0c1c556e560bf76201))
* update seqkit 0.14.0 --> 2.1.0 ([cb99be4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/cb99be4361fbe232be214f0c1c556e560bf76201))
* update snakemake to version 6.13.1 in main environment ([43251f0](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/43251f09945d78cd3a7952d24256e8c703a915f9))
* update tqdm 4.59 --> 4.62 ([485c94e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/485c94ed732057a572f6a3e631aa5618f286e9cf))


### Documentation

* add extra primer keywords ([422f192](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/422f19238d4f91155bd9e6aafba3a8c1dd56aa6a))
* fix typo in docs and update download instructions ([9f239e3](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/9f239e3bf778f2df103a2fc51b0343cd68261b8d))
* update readme with correct primer-keywords ([17878f8](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/17878f84d18774e3e0408dcf13bc555b45b3171c))
* update snakemake version label ([f1a474d](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/f1a474df2cb7556768aabd839a2d9e22910719cd))

### [0.4.2](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.4.1...v0.4.2) (2021-12-08)


### Dependencies

* remove manual installation of pango-designations repo ([5527d53](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/5527d53f554203c5ac59d45df124247ab6eb5396))
* update AmpliGone to version 0.3.3 ([b97a32e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/b97a32ea64054d1581761a9c0e82d0f31b89039e))
* update nextclade typingtool to version 1.6.0 ([303d4d9](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/303d4d9ebc7a6fff966ef4f7a4bd1c1a130c84f8))

### [0.4.1](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.4.0...v0.4.1) (2021-11-26)


### Dependencies

* update AmpliGone to version 0.3.2 ([d47f62b](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/d47f62bd7c8eb632e1f492492efcc76a95ff0ea4))

## [0.4.0](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.3.0...v0.4.0) (2021-11-25)


### Features

* allow SARS2seq to update itself to latest release ([41b711b](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/41b711b6de543d5a6bae2ec6aa10a8afb2ad4917))
* override consensus-index with overlap coordinates index ([b81fd1c](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/b81fd1cf16e95bc113f5598b3094b0262cce1b95))


### Bug Fixes

* Remove --end-bonus and -k28, update -A4 ([e7f8bf7](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/e7f8bf7000f223761df736a83f2bf4b18c81d60c))
* update sars-cov-2 genomic features GFF to include missing ORF ([4399cf4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/4399cf4ba6ef45a68225dad2c1139b74348bf1cc))


### Performance Improvements

* Update minimap2 settings for nanopore data ([75bb475](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/75bb47546a4a108a003c7e694df9ba58972a80bc))


### Documentation

* add (auto-)updating behaviour and "skip-updates" flag to documentation ([43c3655](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/43c36558dfd9f7dd1cf910005ae91facd8c2c439))


### Dependencies

* pin TrueConsense to version 0.3.0 ([b6dfd46](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/b6dfd46bf8ad56b87dca6a528efec51d733172a3))
* update AmpliGone to version 0.3.0 ([f028403](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/f0284030361b5d0bcf0153ea02715efe758c5b35))
* update AmpliGone to version 0.3.1 ([29f1516](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/29f15168fd21c69b7ccdc8c1917da53f10da0253))


### Styles

* remove trailing space ([9933ab4](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/9933ab4c7d27535f6688bc39fc164ebc0eade6d2))

## [0.3.0](https://www.github.com/RIVM-bioinformatics/SARS2seq/compare/v0.2.2...v0.3.0) (2021-11-09)


### Features

* flexible memory requirements per workflow job ([dd1c89e](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/dd1c89e2ddf6f3b7bec0ca48452aabce69da2942))


### Bug Fixes

* use 'map-ont' mm2 preset in nanopore workflow ([87e056c](https://www.github.com/RIVM-bioinformatics/SARS2seq/commit/87e056c03329638ec95e00533e2e16d80424ae52))
