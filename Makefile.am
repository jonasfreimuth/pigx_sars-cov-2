ACLOCAL_AMFLAGS = -I pigx-common/common/m4

nodist_bin_SCRIPTS = \
  pigx-sars-cov-2

dist_pkglibexec_SCRIPTS = \
  snakefile.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS = \
  scripts/create_summary_table.R \
  scripts/create_sample_quality_table.R \
  scripts/deconvolution.R \
  scripts/deconvolution_funs.R \
  scripts/download_databases.sh \
  scripts/fun_index.R \
  scripts/generateNavigation.R \
  scripts/mutation_regression.R \
  scripts/overview_QC_table.R \
  scripts/parse_vep.py \
  scripts/pooling.R \
  scripts/pred_mutation_increase.R \
  scripts/renderReport.R \
  scripts/table_extraction.R \
  scripts/vcfTocsv.py

pkglibexec_snakefilescriptsdir = $(pkglibexecdir)/snakefile_scripts
dist_pkglibexec_snakefilescripts_SCRIPTS = \
  snakefile_scripts/rule_render_qc_report.py \
  snakefile_scripts/download_tarball.py \
  snakefile_scripts/rule_download_kraken_db.py \
  snakefile_scripts/rule_download_krona_db.py \
  snakefile_scripts/rule_download_vep_db.py

pkglibexec_reportsdir = $(pkglibexecdir)/scripts/report_scripts/
dist_pkglibexec_reports_SCRIPTS = \
  scripts/report_scripts/index.Rmd \
  scripts/report_scripts/qc_report_per_sample.Rmd \
  scripts/report_scripts/taxonomic_classification.Rmd \
  scripts/report_scripts/variantreport_per_sample.Rmd

dist_pkgdata_DATA =            \
  etc/pretty.txt \
  etc/settings.yaml \
  etc/sample_sheet.csv.example \
  images/Logo_PiGx.png

dist_doc_DATA = 									\
  README.md

EXTRA_DIST = \
  VERSION \
  META.in \
  $(wildcard tests/vep_test/**) \
  $(wildcard tests/full/**) \
  $(wildcard tests/start_bam/**)


AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1

TESTS = \
  tests/vep_test/vep.sh \
  tests/full/test.sh \
  tests/start_bam/test.sh

clean-local:
	-rm -rf .snakemake
	-rm -rf tests/**/input
	-rm -rf tests/**/output
	-rm -rf tests/**/config.json
	-rm -rf tests/**/*.log

CLEANFILES = $(nodist_bin_SCRIPTS) config.json

integration:
	$(AM_TESTS_ENVIRONMENT) $(abs_top_srcdir)/pigx-sars-cov-2 --printshellcmds \
		-s tests/full/settings.yaml \
		tests/full/sample_sheet.csv

# Build docker image with Guix
docker: dist
	guix pack -C none -e '(load "guix.scm")' --with-source=pigx_sars-cov-2-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils bash
