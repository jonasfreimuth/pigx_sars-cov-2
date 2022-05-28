# PiGx SARS-CoV-2 wastewater sequencing pipeline
#
# Copyright © 2021 Akalin lab.
#
# This file is part of the PiGx SARS-CoV-2 wastewater sequencing pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Snakefile for PiGx SARS-CoV-2 wastewater sequencing pipeline
"""

import os
import csv
import yaml
from itertools import chain
import re
import inspect
from pathlib import Path

SAMPLE_SHEET_CSV = config["locations"]["sample-sheet"]
MUTATION_SHEET_CSV = config["locations"]["mutation-sheet"]
READS_DIR = config["locations"]["reads-dir"]
REFERENCE_FASTA = config["locations"]["reference-fasta"]
AMPLICONS_BED = config["locations"]["amplicons-bed"]
MUTATIONS_BED = config["locations"]["mutations-bed"]
VEP_DB = config["locations"]["vep-db-dir"]
OUTPUT_DIR = config["locations"]["output-dir"]

# TODO: get default read length from multiqc
READ_LENGTH = config["trimming"]["read-length"]
CUT_OFF = config["trimming"]["cut-off"]

MUTATION_COVERAGE_THRESHOLD = config["reporting"]["mutation-coverage-threshold"]

INDEX_DIR = os.path.join(OUTPUT_DIR, "index")
TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, "trimmed_reads")
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")
MAPPED_READS_DIR = os.path.join(OUTPUT_DIR, "mapped_reads")
VARIANTS_DIR = os.path.join(OUTPUT_DIR, "variants")
MUTATIONS_DIR = os.path.join(OUTPUT_DIR, "mutations")
COVERAGE_DIR = os.path.join(OUTPUT_DIR, "coverage")
REPORT_DIR = os.path.join(OUTPUT_DIR, "report")
FASTQC_DIR = os.path.join(REPORT_DIR, "fastqc")
MULTIQC_DIR = os.path.join(REPORT_DIR, "multiqc")
SCRIPTS_DIR = os.path.join(config["locations"]["pkglibexecdir"], "scripts/")
TMP_DIR = os.path.join(config["locations"]["output-dir"], "pigx_work")

if os.getenv("PIGX_UNINSTALLED"):
    LOGO = os.path.join(config["locations"]["pkgdatadir"], "images/Logo_PiGx.png")
else:
    LOGO = os.path.join(config["locations"]["pkgdatadir"], "Logo_PiGx.png")


def toolArgs(name):
    if "args" in config["tools"][name]:
        return config["tools"][name]["args"]
    else:
        return ""


def tool(name):
    cmd = config["tools"][name]["executable"]
    return cmd + " " + toolArgs(name)


BWA_EXEC = tool("bwa")
FASTP_EXEC = tool("fastp")
FASTQC_EXEC = tool("fastqc")
GUNZIP_EXEC = tool("gunzip")
GZIP_EXEC = tool("gzip")
MULTIQC_EXEC = tool("multiqc")
IMPORT_TAXONOMY_EXEC = tool("import_taxonomy")
LOFREQ_EXEC = tool("lofreq")
PYTHON_EXEC = tool("python")
RSCRIPT_EXEC = tool("Rscript")
SAMTOOLS_EXEC = tool("samtools")
VEP_EXEC = tool("vep")
IVAR_EXEC = tool("ivar")

## Load sample sheet
with open(SAMPLE_SHEET_CSV, "r") as fp:
    rows = [row for row in csv.reader(fp, delimiter=",")]
    header = rows[0]
    rows = rows[1:]
    SAMPLE_SHEET = [dict(zip(header, row)) for row in rows]

# Convenience function to access fields of sample sheet columns that
# match the predicate.  The predicate may be a string.
def lookup(column, predicate, fields=[]):
    if inspect.isfunction(predicate):
        records = [line for line in SAMPLE_SHEET if predicate(line[column])]
    else:
        records = [line for line in SAMPLE_SHEET if line[column] == predicate]
    return [record[field] for record in records for field in fields]


SAMPLES = [line["name"] for line in SAMPLE_SHEET]

targets = {
    "help": {"description": "Print all rules and their descriptions.", "files": []},
    "final_reports": {
        "description": "Produce a comprehensive report. This is the default target.",
        "files": (
            expand(
                os.path.join(REPORT_DIR, "{sample}.qc_report_per_sample.html"),
                sample=SAMPLES,
            )
            + expand(
                os.path.join(REPORT_DIR, "{sample}.variantreport_p_sample.html"),
                sample=SAMPLES,
            )
            + [os.path.join(REPORT_DIR, "index.html")]
        ),
    },
    "lofreq": {
        "description": "Call variants and produce .vcf file and overview .csv file.",
        "files": (
            expand(os.path.join(VARIANTS_DIR, "{sample}_snv.csv"), sample=SAMPLES)
        ),
    },
    "multiqc": {
        "description": "Create MultiQC reports for including raw and trimmed reads.",
        "files": (
            expand(
                os.path.join(MULTIQC_DIR, "{sample}", "multiqc_report.html"),
                sample=SAMPLES,
            )
        ),
    },
}
selected_targets = ["final_reports"]
OUTPUT_FILES = list(
    chain.from_iterable([targets[name]["files"] for name in selected_targets])
)


rule all:
    input:
        OUTPUT_FILES,


# Record any existing output files, so that we can detect if they have
# changed.
expected_files = {}


onstart:
    if OUTPUT_FILES:
        for name in OUTPUT_FILES:
            if os.path.exists(name):
                expected_files[name] = os.path.getmtime(name)


# Print generated target files.
onsuccess:
    if OUTPUT_FILES:
        # check if any existing files have been modified
        generated = []
        for name in OUTPUT_FILES:
            if (
                name not in expected_files
                or os.path.getmtime(name) != expected_files[name]
            ):
                generated.append(name)
        if generated:
            print("The following files have been generated:")
            for name in generated:
                print("  - {}".format(name))

# Trimming in three steps: general by qual and cutoff, get remaining adapters out, get remaining primers out


rule bwa_index:
    input:
        REFERENCE_FASTA,
    output:
        ref=os.path.join(INDEX_DIR, os.path.basename(REFERENCE_FASTA)),
        index=os.path.join(
            INDEX_DIR, "{}.bwt".format(os.path.basename(REFERENCE_FASTA))
        ),
    log:
        os.path.join(LOG_DIR, "bwa_index.log"),
    shell:
        """
        mkdir -p {INDEX_DIR};
        ln -sf {input} {INDEX_DIR};
        cd {INDEX_DIR};
        {BWA_EXEC} index {output.ref} >> {log} 2>&1
        """


rule samtools_sort_preprimertrim:
    input:
        os.path.join(READS_DIR, "{sample}.bam"),
    output:
        os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bam"),
    log:
        os.path.join(LOG_DIR, "samtools_sort_{sample}.log"),
    shell:
        "{SAMTOOLS_EXEC} sort -o {output} {input} >> {log} 2>&1"


rule samtools_index_preprimertrim:
    input:
        os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bam"),
    output:
        os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bai"),
    log:
        os.path.join(LOG_DIR, "samtools_index_{sample}.log"),
    shell:
        "{SAMTOOLS_EXEC} index {input} {output} >> {log} 2>&1"


# WIP create a dummy entry if no variant is found - use this as long as the input-function solution doesn't work
def no_variant_vep(sample, lofreq_output):
    content = open(lofreq_output.format(sample=sample), "r").read()
    if re.findall("^NC", content, re.MULTILINE):  # regex ok or not?
        # trigger vep path
        print("File can be used for downstream processing")
    else:
        # write smth so that vep does not crash - deal with everything later in the variant_report
        print("adding dummy entry to vcf file, because no variants were found")
        open(lofreq_output.format(sample=sample), "a").write(
            "NC_000000.0\t00\t.\tA\tA\t00\tPASS\tDP=0;AF=0;SB=0;DP4=0,0,0,0"
        )


# TODO it should be possible to add customized parameter
rule lofreq:
    input:
        aligned_bam=os.path.join(
            MAPPED_READS_DIR, "{sample}_aligned_sorted.bam"
        ),
        aligned_bai=os.path.join(
            MAPPED_READS_DIR, "{sample}_aligned_sorted.bai"
        ),
        ref=os.path.join(INDEX_DIR, "{}".format(os.path.basename(REFERENCE_FASTA))),
    output:
        vcf=os.path.join(VARIANTS_DIR, "{sample}_snv.vcf"),
    log:
        os.path.join(LOG_DIR, "lofreq_{sample}.log"),
    run:
        shell(
            "{LOFREQ_EXEC} call -f {input.ref} -o {output} --verbose {input.aligned_bam} >> {log} 2>&1"
        )
        # WIP create a dummy entry if no variant is found - use this as long as the input-function solution doesn't work
        no_variant_vep(wildcards.sample, output.vcf)


rule vcf2csv:
    input:
        os.path.join(VARIANTS_DIR, "{sample}_snv.vcf"),
    output:
        os.path.join(VARIANTS_DIR, "{sample}_snv.csv"),
    params:
        script=os.path.join(SCRIPTS_DIR, "vcfTocsv.py"),
    log:
        os.path.join(LOG_DIR, "vcf2csv_{sample}.log"),
    shell:
        "{PYTHON_EXEC} {params.script} {input} >> {log} 2>&1"


rule vep:
    input:
        os.path.join(VARIANTS_DIR, "{sample}_snv.vcf"),
    output:
        os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2.txt"),
    params:
        species="sars_cov_2",
    log:
        os.path.join(LOG_DIR, "vep_{sample}.log"),
    shell:
        """
        {VEP_EXEC} --verbose --offline --dir_cache {VEP_DB} --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing\
        --distance 5000 --mane --protein --species {params.species} --symbol --transcript_version --tsl\
        --input_file {input} --output_file {output} >> {log} 2>&1
        """


rule parse_vep:
    input:
        os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2.txt"),
    output:
        os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2_parsed.txt"),
    params:
        script=os.path.join(SCRIPTS_DIR, "parse_vep.py"),
    log:
        os.path.join(LOG_DIR, "parse_vep_{sample}.log"),
    shell:
        "{PYTHON_EXEC} {params.script} {input} {output} >> {log} 2>&1"


# TODO: change amplicon naming to mutation site since it is misleading
rule samtools_bedcov:
    input:
        mutations_bed=MUTATIONS_BED,
        aligned_bam=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bam"),
        aligned_bai=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bai"),
    output:
        os.path.join(COVERAGE_DIR, "{sample}_amplicons.csv"),
    log:
        os.path.join(LOG_DIR, "samtools_bedcov_{sample}.log"),
    shell:
        "{SAMTOOLS_EXEC} bedcov {input.mutations_bed} {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule samtools_coverage:
    input:
        aligned_bam=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bam"),
        aligned_bai=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bai"),
    output:
        os.path.join(COVERAGE_DIR, "{sample}_coverage.csv"),
    log:
        os.path.join(LOG_DIR, "samtools_coverage_{sample}.log"),
    shell:
        "{SAMTOOLS_EXEC} coverage {input.aligned_bam} > {output} 2>> {log} 3>&2"


# TODO: change amplicon naming to mutation site because it is misleading
rule get_qc_table:
    input:
        coverage_csv=os.path.join(COVERAGE_DIR, "{sample}_coverage.csv"),
        amplicon_csv=os.path.join(COVERAGE_DIR, "{sample}_amplicons.csv"),
    output:
        os.path.join(COVERAGE_DIR, "{sample}_merged_covs.csv"),
    params:
        script=os.path.join(SCRIPTS_DIR, "get_qc_table.py"),
    log:
        os.path.join(LOG_DIR, "get_qc_table_{sample}.log"),
    shell:
        "{PYTHON_EXEC} {params.script} {input.coverage_csv} {input.amplicon_csv} {output} >> {log} 2>&1"


rule generate_navbar:
    input:
        script=os.path.join(SCRIPTS_DIR, "generateNavigation.R"),
    output:
        os.path.join(REPORT_DIR, "_navbar.html"),
    params:
        report_scripts_dir=os.path.join(SCRIPTS_DIR, "report_scripts"),
    log:
        os.path.join(LOG_DIR, "generate_navigation.log"),
    shell:
        "{RSCRIPT_EXEC} {input.script} \
        {params.report_scripts_dir} {SAMPLE_SHEET_CSV} {output} > {log} 2>&1"


rule render_variant_report:
    input:
        vep=os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2_parsed.txt"),
        snv=os.path.join(VARIANTS_DIR, "{sample}_snv.csv"),
        deconvolution_functions=os.path.join(SCRIPTS_DIR, "deconvolution.R"),
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "variantreport_p_sample.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
    output:
        varreport=os.path.join(REPORT_DIR, "{sample}.variantreport_p_sample.html"),
        mutations=os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"),
        variants=os.path.join(VARIANTS_DIR, "{sample}_variants.csv"),
    log:
        os.path.join(LOG_DIR, "reports", "{sample}_variant_report.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
        {input.report} {output.varreport} {input.header} \
        '{{\
          "sample_name":  "{wildcards.sample}",  \
          "output_dir": "{OUTPUT_DIR}",      \
          "vep_file":     "{input.vep}",         \
          "snv_file":     "{input.snv}",         \
          "sample_sheet": "{SAMPLE_SHEET_CSV}",  \
          "mutation_sheet": "{MUTATION_SHEET_CSV}", \
          "deconvolution_functions": "{input.deconvolution_functions}", \
          "logo": "{LOGO}" \
        }}' > {log} 2>&1
        """


rule render_qc_report:
    input:
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "qc_report_per_sample.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
        coverage=os.path.join(COVERAGE_DIR, "{sample}_merged_covs.csv")
    output:
        os.path.join(REPORT_DIR, "{sample}.qc_report_per_sample.html"),
    params:
        multiqc_rel_path=lambda wildcards, input: input.multiqc[len(REPORT_DIR) + 1 :],
    log:
        os.path.join(LOG_DIR, "reports", "{sample}_qc_report.log"),
    shell:
        """{RSCRIPT_EXEC} {input.script} \
        {input.report} {output} {input.header} \
        '{{\
          "sample_name": "{wildcards.sample}",  \
          "coverage_file": "{input.coverage}",   \
          "logo": "{LOGO}" \
        }}' > {log} 2>&1"""


rule create_variants_summary:
    input:
        script=os.path.join(SCRIPTS_DIR, "creating_variants_summary_table.R"),
        files=expand(
            os.path.join(VARIANTS_DIR, "{sample}_variants.csv"), sample=SAMPLES
        ),
    output:
        os.path.join(VARIANTS_DIR, "data_variant_plot.csv"),
    log:
        os.path.join(LOG_DIR, "create_variants_summary.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} "{VARIANTS_DIR}" {output} > {log} 2>&1
        """


rule create_mutations_summary:
    input:
        script=os.path.join(SCRIPTS_DIR, "creating_mutation_summary_table.R"),
        files=expand(
            os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"), sample=SAMPLES
        ),
    output:
        os.path.join(VARIANTS_DIR, "data_mutation_plot.csv"),
    log:
        os.path.join(LOG_DIR, "create_mutations_summary.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} "{MUTATIONS_DIR}" {output} > {log} 2>&1
        """


# TODO integrate the output of fastp.json to get the number of raw and trimmed reads
rule create_overviewQC_table:
    input:
        script=os.path.join(SCRIPTS_DIR, "overview_QC_table.R"),
        cov_summary=expand(
            os.path.join(COVERAGE_DIR, "{sample}_merged_covs.csv"), sample=SAMPLES
        ),
    output:
        os.path.join(OUTPUT_DIR, "overview_QC.csv"),
    log:
        os.path.join(LOG_DIR, "create_overviewQC_table.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} {SAMPLE_SHEET_CSV} {output} {READS_DIR} {TRIMMED_READS_DIR} {MAPPED_READS_DIR} {COVERAGE_DIR} > {log} 2>&1
        """


rule render_index:
    input:
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "index.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
        variants=os.path.join(VARIANTS_DIR, "data_variant_plot.csv"),
        mutations=os.path.join(VARIANTS_DIR, "data_mutation_plot.csv"),
        overviewQC=os.path.join(OUTPUT_DIR, "overview_QC.csv"),
    params:
        fun_cvrg_scr=os.path.join(SCRIPTS_DIR, "sample_coverage_score.R"),
        fun_lm=os.path.join(SCRIPTS_DIR, "pred_mutation_increase.R"),
        fun_tbls=os.path.join(SCRIPTS_DIR, "table_extraction.R"),
        fun_pool=os.path.join(SCRIPTS_DIR, "pooling.R"),
    output:
        report=os.path.join(REPORT_DIR, "index.html"),
        tbl_mut_count=os.path.join(OUTPUT_DIR, "mutations_counts.csv"),
        tbl_lm_res=os.path.join(OUTPUT_DIR, "linear_regression_results.csv"),
    log:
        os.path.join(LOG_DIR, "reports", "index.log"),
    shell:
        """{RSCRIPT_EXEC} {input.script} \
        {input.report} {output.report} {input.header}   \
        '{{                                      \
          "variants_csv": "{input.variants}",   \
          "mutations_csv": "{input.mutations}", \
          "coverage_dir": "{COVERAGE_DIR}",\
          "sample_sheet": "{SAMPLE_SHEET_CSV}",  \
          "mutation_sheet": "{MUTATION_SHEET_CSV}", \
          "mutation_coverage_threshold": "{MUTATION_COVERAGE_THRESHOLD}", \
          "logo": "{LOGO}", \
          "fun_cvrg_scr": "{params.fun_cvrg_scr}", \
          "fun_lm": "{params.fun_lm}", \
          "fun_tbls": "{params.fun_tbls}", \
          "fun_pool": "{params.fun_pool}", \
          "overviewQC": "{input.overviewQC}", \
          "output_dir": "{OUTPUT_DIR}" \
        }}' > {log} 2>&1"""
