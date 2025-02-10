## Bump new analysis version

Read the docs  
https://gitlab.curie.fr/dm-toolbox/deploy-geniac/-/blob/devel/docs/gitlab-ci-data-analysis.md

### Code and versioning (on the `devel` branch)

- [ ] The version number (MAJOR.UPDATE.MINOR) is updated (`ǹextflow.config`) and the `dev` suffix is removed 
- [ ] The CHANGELOG is updated and summarize the issues related to each modification
- [ ] The README page is up to date
- [ ] The `parameters.settings.json`file is updated with the last parameters
- [ ] The `usage.md` documentation file is updated
- [ ] The `nf-modules` folder is updated and compared to the modules repository
- [ ] The `recipes`folder is up-to-date with the last conda environment files
- [ ] The annotation config file is up-to-date
- [ ] The milestone related to this new version is completed
- [ ] The .gitlab-ci.yml is the latest version
- [ ] The output files are exported and organized as expected
- [ ] The MultiQC report is validated
- [ ] Check if the code can be 100% resume

### Lint

- [ ] The geniac linter is running without issue
- [ ] The CI is running

### Testing

- [ ] The test profile is running with conda
- [ ] The test profile is running with singularity

### Bump

- [ ] Push the new version on the master/main branch
- [ ] Create a new tag/release from the master branch for this version
- [ ] Update the `devel`branch to be ready for the next version and add again the `dev`flag