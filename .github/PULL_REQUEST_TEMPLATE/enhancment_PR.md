---
name: Enhancement PR
about: Submit a PR for an enhancement
title: 'ENH: '
labels: enhancement
assignees: ''

---
### About this repo
- `q2-moshpit` is a `qiime2` plugin that offers various analyses for shotgun metagenome data.
- <!--- What part of the repo is the current PR concerned with? E.g. The current PR addressed a bug in action X, which runs a Y analysis on the inputs. -->

### What's new
- <!--- Describe what was changed in the code and why it is useful or necessary-->
- <!--- Does this PR fully address an existing issue? If so write: Closes #issue_number -->

### Additional Info
<!--- Is the PR blocked by another PR? If so, disclose it here. You can use the syntax user/repo_name/pull/PR_number to reference PRs in other repos. To reference PRs in the same repo simply use #PR_number. Do so inside the ```[tasklist]``` context as shown below.
```[tasklist]
### Blocked by
- [ ] user/repo_name/pull/PR_number
- [ ] #PR_number
```
If you also what the CI tool to block merging of the PR before another you can use the following:
- blocked by #PR_number
- merge after user/repo_name#PR_number
- dependent on #PR_number
-->

### Set up an environment
<!---
- The following commands should get the reviewer a working environment where they can test the PR changes. 
- Keep in mind that if you PR depends on newer versions of the dependencies you will have to install these manually, e.g: pip install git+https://github.com/username/repository.git or pip install -e . your_local_dependency
-->
```bash
# For linux: 
# export MY_OS="linux"
# For mac:
export MY_OS="osx" 
wget "https://data.qiime2.org/distro/shotgun/qiime2-shotgun-2023.9-py38-"$MY_OS"-conda.yml"
conda env create -n q2-shotgun --file qiime2-shotgun-2023.9-py38-osx-conda.yml
rm "qiime2-shotgun-2023.9-py38-"$MY_OS"-conda.yml"
```

### Run it locally 
1. Clone the repo and checkout the PR branch (skip this step if you already have a local copy of the code):
```bash
# Remove q2-moshpit so you can install your local version.
conda activate q2-shotgun
conda remove q2-moshpit

# clone from GitHub
git clone git@github.com:bokulich-lab/q2-moshpit.git
cd q2-moshpit
gh pr checkout <PR_number>
pip install -e .
```
<!---
- The PR_number will be created after you submit the PR, therefore it can only be set after, by editing the PR message.
- Make sure you have gh (GitHub CLI) installed in your environment.
-->

2. Let's get you some data to play with: 
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
<!---Feel free to eliminate sections that are not relevant to your PR.-->
- [ ] Fill in *About this repo* section.
- [ ] Fill in *Whats new* section. Erase any bullet points that are not used.
- [ ] In the *Set up an environment* Make sure to adjust the environment to the latest version possible of all involved dependencies.
- [ ] In *Run it locally*, **step 1**, make sure to write the PR number after you have submitted the PR.
- [ ] In *Run it locally*, **step 2**, fill in the code to fetch the inputs to test the code.
- [ ] In *Run it locally*, **step 3**, fill in the code to test the changes. 
