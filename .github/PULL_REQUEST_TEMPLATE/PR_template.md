### What's new
<!---
- Describe what was changed in the code and why it is useful or necessary
- Does this PR fully address an existing issue? If so write: Closes #issue_number 
-->

### Additional Info
<!--- Is the PR blocked by another PR? If so, disclose it here. You can use the syntax user/repo_name/pull/PR_number to reference PRs in other repos. To reference PRs in the same repo simply use #PR_number. Do so inside the ```[tasklist]``` context as shown below.

```[tasklist]
### Blocked by
- [ ] user/repo_name#PR_number
- [ ] #PR_number
```
-->

<!---
If you also what the CI tool to block merging of the PR before another use the following syntax:
- blocked by #PR_number
- merge after user/repo_name#PR_number
- dependent on #PR_number

Feel free to leave this section as a comment. The CI will still pick it up. 
-->

### Run it locally 
1. Assuming you already a working virtual environment.
```bash
# Remove q2-moshpit so you can install your local version.
conda activate <you_env_name>
conda remove q2-moshpit
```

```bash
# clone from GitHub
cd <clone_here>
git clone git@github.com:bokulich-lab/q2-moshpit.git
cd q2-moshpit
pip install -e .
```

```bash
# checkout the PR branch
PR=<PR_number>
git fetch origin pull/${PR}/head:pr-${PR}
git checkout pr-${PR}
```
<!---
- The PR_number will be created after you submit the PR, therefore it can only be set after, by editing the PR message.
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

### TODO
<!---Feel free to eliminate sections that are not relevant to your PR.-->
- [ ] Fill in *Whats new* section. Erase any bullet points that are not used.
- [ ] In *Run it locally*, **step 1**, make sure to write the PR number after you have submitted the PR.
- [ ] In *Run it locally*, **step 2**, fill in the code to fetch the inputs to test the code.
- [ ] In *Run it locally*, **step 3**, fill in the code to test the changes. 
