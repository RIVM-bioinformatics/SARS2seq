# Changelog

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
