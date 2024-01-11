---
name: PR template
about: Submit a PR
title: 'BUG: '
labels: bug
assignees: ''

---
### About this repo
- `q2-moshpit` is a `qiime2` plugin that offers various analyses for shotgun metagenome data.
- <!--- What part of the repo is the current PR concerned with? E.g. The current PR addressed a bug in action X, which runs a Y analysis on the inputs. -->

### What's new
- <!--- Does this PR fully address an existing issue? If so write: Closes #issue_number -->
- <!--- Describe what was changed in the code -->
- <!--- Describe why the changes are useful or necessary -->
- <!--- Is the PR blocked by another PR? If so, disclose it here. You can use the syntax user/repo_name/pull/PR_number to reference PRs in other repos. To reference PRs in the same repo simply use #PR_number -->

### Set up an environment
<!--- The following commands should get the reviewer a working environment where they can test the PR changes. -->
```bash
mamba create -yn q2-shotgun \
    -c conda-forge -c bioconda -c https://packages.qiime2.org/qiime2/2023.5/tested -c defaults \
    q2cli q2-moshpit gh

conda run -n q2-shotgun \
    pip install --no-deps --force-reinstall git+https://github.com/misialq/quast.git@issue-230
```

### Run it locally 
1. Clone the repo and checkout the PR branch (skip this step if you already have a local copy of the PR branch):
```bash
git clone git@github.com:bokulich-lab/q2-moshpit.git
gh pr checkout <PR_number>
```
<!---
- The PR_number will be created after you submit the PR, therefore it can only be set after, by editing the PR message.
- Make sure you have gh (GitHub CLI) installed in your environment.
-->

2. Lets get you some data to play with: 
<!---In the next steps provide terminal commands that get the reviewer the necessary inputs to run a working example.-->
```bash
<your code here>
```
<!---
Example:
```
wget https://scop.berkeley.edu/downloads/scopeseq-2.07/astral-scopedom-seqres-gd-sel-gs-bib-40-2.07.fa
mkdir sequences
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $0}}' astral-scopedom-seqres-gd-sel-gs-bib-40-2.07.fa > sequences/protein-sequences.fasta
qiime tools import --input-path sequences --output-path sequences.qza --type FeatureData\[ProteinSequence\]
```
-->

3. Test it out!
```bash
<your code here>
```
<!---
Example:
```
qiime moshpit build-diamond-db --i-sequences sequences.qza --o-diamond-db custome_diamond.qza --verbose
```
-->

### Running the tests
```bash
pytest -W ignore -vv --pyargs q2_moshpit
```

### TODO
- [ ] Fill in *About this repo* section.
- [ ] Fill in *Whats new* section. Erase any bullet points that are not used.
- [ ] In the *Set up an environment* Make sure to adjust the environment to the latest version possible of all involved dependencies.
- [ ] In *Run it locally*, **step 1**, make sure to write the PR number after you have submitted the PR.
- [ ] In *Run it locally*, **step 2**, fill in the code to fetch the inputs to test the code.
- [ ] In *Run it locally*, **step 3**, fill in the code to test the changes. 
