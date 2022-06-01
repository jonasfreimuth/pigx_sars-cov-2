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


# Helper functions


def toolArgs(name):
    """
    Helper function to retrieve all preset arguments for a given tool.
    """
    if "args" in config["tools"][name]:
        return config["tools"][name]["args"]
    else:
        return ""


def tool(name):
    """
    Helper function to bundle a tool with its preset arguments.
    """
    cmd = config["tools"][name]["executable"]
    return cmd + " " + toolArgs(name)


# Convenience function to access fields of sample sheet columns that
# match the predicate.  The predicate may be a string.
def lookup(column, predicate, fields=[]):
    """
    Convenience function to access fields of sample sheet columns that
    match the predicate.  The predicate may be a string.
    """
    if inspect.isfunction(predicate):
        records = [line for line in SAMPLE_SHEET if predicate(line[column])]
    else:
        records = [line for line in SAMPLE_SHEET if line[column] == predicate]
    return [record[field] for record in records for field in fields]


def fastq_ext(fastq_file):
    "Function to determine the fastq file extension"
    root, ext = os.path.splitext(fastq_file)
    if ext == ".gz":
        root_root, root_ext = os.path.splitext(root)
        ext = "".join([root_ext, ext])
    return ext


# Input functions


# function to pass read files to trim/filter/qc improvement
def trim_reads_input(args):
    """
    Get a list of all files related to a sample from the sample sheet. Helps in
    working with both single and paired end data.
    """
    sample = args[0]
    return [
        os.path.join(READS_DIR, f)
        for f in lookup("name", sample, ["reads", "reads2"])
        if f
    ]


def map_input(args):
    """
    Function to return the trimmed read files belonging to a sample, independent
    of their end mode (single vs paired).
    """
    sample = args[0]
    reads_files = [
        os.path.join(READS_DIR, f)
        for f in lookup("name", sample, ["reads", "reads2"])
        if f
    ]
    if len(reads_files) > 1:
        return [
            os.path.join(
                TRIMMED_READS_DIR, "{sample}_trimmed_R1.fastq.gz".format(sample=sample)
            ),
            os.path.join(
                TRIMMED_READS_DIR, "{sample}_trimmed_R2.fastq.gz".format(sample=sample)
            ),
        ]
    elif len(reads_files) == 1:
        return [
            os.path.join(
                TRIMMED_READS_DIR, "{sample}_trimmed.fastq.gz".format(sample=sample)
            )
        ]


# dynamically define the multiqc input files created by FastQC and fastp
# TODO add kraken reports per sample
def multiqc_input(args):
    """
    Dynamically define the multiqc input files created by FastQC and fastp for
    each sample.
    """
    sample = args[0]
    reads_files = [
        os.path.join(READS_DIR, f)
        for f in lookup("name", sample, ["reads", "reads2"])
        if f
    ]
    # read_num is either ["_R1", "_R2"] or [""] depending on number of read files
    read_num = [
        "_R" + str(f) if len(reads_files) > 1 else ""
        for f in range(1, len(reads_files) + 1)
    ]
    se_or_pe = ["pe" if len(reads_files) > 1 else "se"]
    files = [
        # fastp on raw files
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}_{end}_fastp.html"),
            sample=sample,
            end=se_or_pe,
        ),
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}_{end}_fastp.json"),
            sample=sample,
            end=se_or_pe,
        ),
        # fastqc on raw files
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}{read_num}_fastqc.html"),
            sample=sample,
            read_num=read_num,
        ),
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}{read_num}_fastqc.zip"),
            sample=sample,
            read_num=read_num,
        ),
        # fastqc after trimming
        expand(
            os.path.join(
                FASTQC_DIR, "{sample}", "{sample}_trimmed{read_num}_fastqc.html"
            ),
            sample=sample,
            read_num=read_num,
        ),
        expand(
            os.path.join(
                FASTQC_DIR, "{sample}", "{sample}_trimmed{read_num}_fastqc.zip"
            ),
            sample=sample,
            read_num=read_num,
        ),
        # fastqc after primer trimming
        expand(
            os.path.join(
                FASTQC_DIR,
                "{sample}",
                "{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.html",
            ),
            sample=sample,
        ),
        expand(
            os.path.join(
                FASTQC_DIR,
                "{sample}",
                "{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.zip",
            ),
            sample=sample,
        ),
    ]
    return list(chain.from_iterable(files))


# WIP - until then use hack that create a single line in the lofreq output
def vep_input(args):
    """
    Fucntion to skip all VEP rules for a sample in case no variants are found by
    lofreq through creation of dummy VEP output files (?)
    """
    sample = args[0]
    lofreq_output = (
        rules.lofreq.output.vcf
    )  # this requires the file to be there already - I have no idea how to make the decision about the further input when it requires a rule to run beforhand
    logger.info(lofreq_output.format(sample=sample))
    with open(lofreq_output.format(sample=sample), "r") as vcf:
        content = vcf.read()
        if re.findall("^NC", content, re.MULTILINE):  # regex ok or not?
            # trigger vep path
            return [
                os.path.join(
                    VARIANTS_DIR,
                    "{sample}_vep_sarscov2_parsed.txt".format(sample=sample),
                ),
                os.path.join(VARIANTS_DIR, "{sample}_snv.csv".format(sample=sample)),
            ]
        else:
            # skipp execution of all vep related rules and directly have smth that the report can work with
            empty_vep_txt = os.path.join(
                VARIANTS_DIR, "{sample}_vep_sarscov2_empty.txt".format(sample=sample)
            )
            Path(empty_vep_txt).touch()
            empty_snv_csv = os.path.join(
                VARIANTS_DIR, "{sample}_snv_empty.csv".format(sample=sample)
            )
            Path(empty_snv_csv).touch()
            return [empty_vep_txt, empty_snv_csv]


# WIP create a dummy entry if no variant is found - use this as long as the input-function solution doesn't work
def no_variant_vep(sample, lofreq_output):
    """
    Work-around to create dummy entries in lofreq output to ensure smooth VEP
    running in case no variants are found by lofreq (?)
    """
    content = open(lofreq_output.format(sample=sample), "r").read()
    if re.findall("^NC", content, re.MULTILINE):  # regex ok or not?
        # trigger vep path
        logger.info("File can be used for downstream processing")
    else:
        # write smth so that vep does not crash - deal with everything later in the variant_report
        logger.info("adding dummy entry to vcf file, because no variants were found")
        open(lofreq_output.format(sample=sample), "a").write(
            "NC_000000.0\t00\t.\tA\tA\t00\tPASS\tDP=0;AF=0;SB=0;DP4=0,0,0,0")


SAMPLE_SHEET_CSV    = config["locations"]["sample-sheet"]
MUTATION_SHEET_CSV  = config["locations"]["mutation-sheet"]
READS_DIR           = config["locations"]["reads-dir"]
REFERENCE_FASTA     = config["locations"]["reference-fasta"]
AMPLICONS_BED       = config["locations"]["amplicons-bed"]
MUTATIONS_BED       = config["locations"]["mutations-bed"]
KRAKEN_DB           = config["locations"]["kraken-db-dir"]
KRONA_DB            = config["locations"]["krona-db-dir"]
VEP_DB              = config["locations"]["vep-db-dir"]
OUTPUT_DIR          = config["locations"]["output-dir"]

# TODO: get default read length from multiqc
READ_LENGTH      = config['trimming']['read-length']
CUT_OFF          = config['trimming']['cut-off']

MUTATION_DEPTH_THRESHOLD    = config["reporting"]["mutation-depth-threshold"]
MUTATION_COVERAGE_THRESHOLD = config['reporting']['mutation-coverage-threshold']

START_POINT = config["control"]["start"]
TARGETS     = config["control"]["targets"]

RUN_IVAR_PRIMER_TRIMING = config["control"]["run-ivar-primer-triming"]

INDEX_DIR         = os.path.join(OUTPUT_DIR, 'index')
TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'mapped_reads')
VARIANTS_DIR      = os.path.join(OUTPUT_DIR, 'variants')
MUTATIONS_DIR     = os.path.join(OUTPUT_DIR, 'mutations')
KRAKEN_DIR        = os.path.join(OUTPUT_DIR, 'kraken')
COVERAGE_DIR      = os.path.join(OUTPUT_DIR, 'coverage')
REPORT_DIR        = os.path.join(OUTPUT_DIR, 'report')
FASTQC_DIR        = os.path.join(REPORT_DIR, 'fastqc')
MULTIQC_DIR       = os.path.join(REPORT_DIR, 'multiqc')
SCRIPTS_DIR       = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
TMP_DIR           = os.path.join(config['locations']['output-dir'], 'pigx_work')

if os.getenv("PIGX_UNINSTALLED"):
    LOGO = os.path.join(config['locations']['pkgdatadir'], "images/Logo_PiGx.png")
else:
    LOGO = os.path.join(config['locations']['pkgdatadir'], "Logo_PiGx.png")

BWA_EXEC             = tool("bwa")
FASTP_EXEC           = tool("fastp")
FASTQC_EXEC          = tool("fastqc")
GUNZIP_EXEC          = tool("gunzip")
GZIP_EXEC            = tool("gzip")
MULTIQC_EXEC         = tool("multiqc")
IMPORT_TAXONOMY_EXEC = tool("import_taxonomy")
KRAKEN2_EXEC         = tool("kraken2")
LOFREQ_EXEC          = tool("lofreq")
PYTHON_EXEC          = tool("python")
RSCRIPT_EXEC         = tool("Rscript")
SAMTOOLS_EXEC        = tool("samtools")
VEP_EXEC             = tool("vep")
IVAR_EXEC            = tool("ivar")

## Load sample sheet
with open(SAMPLE_SHEET_CSV, 'r') as fp:
  rows =  [row for row in csv.reader(fp, delimiter=',')]
  header = rows[0]; rows = rows[1:]
  SAMPLE_SHEET = [dict(zip(header, row)) for row in rows]

SAMPLES = [line['name'] for line in SAMPLE_SHEET]

targets = {
    'help': {
        'description': "Print all rules and their descriptions.",
        'files': []
    },
    'final_reports': {
        'description': "Produce a comprehensive report. This is the default target.",
        'files': (
            expand(os.path.join(REPORT_DIR, '{sample}.qc_report_per_sample.html'), sample=SAMPLES) +
            expand(os.path.join(REPORT_DIR, '{sample}.variantreport_p_sample.html'), sample=SAMPLES) +
            expand(os.path.join(REPORT_DIR, '{sample}.taxonomic_classification.html'), sample=SAMPLES) +
            expand(os.path.join(REPORT_DIR, '{sample}.Krona_report.html'), sample=SAMPLES) +
            [os.path.join(REPORT_DIR, 'index.html')]
        )
    },
    'lofreq': {
        'description': "Call variants and produce .vcf file and overview .csv file.",
        'files': (
            expand(os.path.join(VARIANTS_DIR, '{sample}_snv.csv'), sample=SAMPLES)
        )
    },
    'multiqc': {
        'description': "Create MultiQC reports for including raw and trimmed reads.",
        'files': (
            expand(os.path.join(MULTIQC_DIR, '{sample}', 'multiqc_report.html'), sample=SAMPLES)
        )
    }
}

selected_targets = config["control"]["targets"]
OUTPUT_FILES = list(chain.from_iterable([targets[name]['files'] for name in selected_targets]))


rule all:
    input: OUTPUT_FILES

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
            if name not in expected_files or os.path.getmtime(name) != expected_files[name]:
                generated.append(name)
        if generated:
            print("The following files have been generated:")
            for name in generated:
                logger.info("  - {}".format(name))

# Trimming in three steps: general by qual and cutoff, get remaining adapters out, get remaining primers out

# TODO the output suffix should be dynamic depending on the input
# TODO with the use of fastp the use of fastqc becomes partly reduntant, fastqc should be removed or adjusted
rule fastp:
    input: trim_reads_input
    output:
        r1 = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R1.fastq.gz"),
        r2 = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R2.fastq.gz"),
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_pe_fastp.html'),
        json = os.path.join(FASTQC_DIR, '{sample}', '{sample}_pe_fastp.json')
    log: os.path.join(LOG_DIR, 'fastp_{sample}.log')
    shell: """
         {FASTP_EXEC} -i {input[0]} -I {input[1]} -o {output.r1} -O {output.r2} --html {output.html} --json {output.json} >> {log}t 2>&1
     """

rule fastp_se:
    input: trim_reads_input
    output:
        r = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed.fastq.gz"),
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_se_fastp.html'),
        json = os.path.join(FASTQC_DIR, '{sample}', '{sample}_se_fastp.json')
    log: os.path.join(LOG_DIR, 'fastp_{sample}.log')
    shell: """
        {FASTP_EXEC} -i {input[0]} -o {output.r} --html {output.html} --json {output.json} >> {log}t 2>&1
    """

rule bwa_index:
    input: REFERENCE_FASTA
    output:
      ref=os.path.join(INDEX_DIR, os.path.basename(REFERENCE_FASTA)),
      index=os.path.join(INDEX_DIR, "{}.bwt".format(os.path.basename(REFERENCE_FASTA)))
    log: os.path.join(LOG_DIR, 'bwa_index.log')
    shell: """
        mkdir -p {INDEX_DIR};
        ln -sf {input} {INDEX_DIR};
        cd {INDEX_DIR};
        {BWA_EXEC} index {output.ref} >> {log} 2>&1
        """

# alignment works with both single and paired-end files
rule bwa_align:
    input:
        fastq = map_input,
        ref = os.path.join(INDEX_DIR, "{}".format(os.path.basename(REFERENCE_FASTA))),
        index = os.path.join(INDEX_DIR, "{}.bwt".format(os.path.basename(REFERENCE_FASTA)))
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    params:
        threads = 4
    log: os.path.join(LOG_DIR, 'bwa_align_{sample}.log')
    shell: "{BWA_EXEC} mem -t {params.threads} {input.ref} {input.fastq} > {output} 2>> {log} 3>&2"

# TODO verify that subsequent tools do not require filtering for proper pairs
# NOTE verification of flags can be done with "samtools view -h -f 4 <file.sam|file.bam> | samtools flagstat -
rule samtools_filter_aligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    params:
        # add 'proper-pair' filter (-f 2) if sample is paired-end
        proper_pair = lambda wc: "-f 2" if len(trim_reads_input(wc))>1 else ""
    log: os.path.join(LOG_DIR, 'samtools_filter_aligned_{sample}.log')
    shell: # exclude (F) reads that are not mapped (4) and supplementary (2048)
        "{SAMTOOLS_EXEC} view -bh {params.proper_pair} -F 4 -F 2048 {input} > {output} 2>> {log} 3>&2"

rule samtools_filter_unaligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    log: os.path.join(LOG_DIR, 'samtools_filter_unaligned_{sample}.log')
    shell: # keep (-f) reads that are unmapped (4)
        "{SAMTOOLS_EXEC} view -bh -f 4 {input} > {output} 2>> {log} 3>&2"

rule samtools_sort_preprimertrim:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    log: os.path.join(LOG_DIR, 'samtools_sort_{sample}.log')
    shell: "{SAMTOOLS_EXEC} sort -o {output} {input} >> {log} 2>&1"

rule samtools_index_preprimertrim:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    log: os.path.join(LOG_DIR, 'samtools_index_{sample}.log')
    shell: "{SAMTOOLS_EXEC} index {input} {output} >> {log} 2>&1"

rule ivar_primer_trim:
    input:
        primers=AMPLICONS_BED,
        aligned_bam=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bam"),
        aligned_bai=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bai"),
    output:
        os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted_primer-trimmed.bam"),
    params:
        output=lambda wildcards, output: os.path.splitext(f"{output}")[0],
    log:
        os.path.join(LOG_DIR, "ivar_{sample}.log"),
    # TODO number parameter should be accessible over settings file
    shell:
        """
        {IVAR_EXEC} trim -b {input.primers} -p {params.output} -i {input.aligned_bam} -q 15 -m 180 -s 4 >> {log} 2>&1
        """


rule samtools_sort_postprimertrim:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bam')
    log: os.path.join(LOG_DIR, 'samtools_sort_{sample}.log')
    shell: "{SAMTOOLS_EXEC} sort -o {output} {input} >> {log} 2>&1"

rule samtools_index_postprimertrim:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bai')
    log: os.path.join(LOG_DIR, 'samtools_index_{sample}.log')
    shell: "{SAMTOOLS_EXEC} index {input} {output} >> {log} 2>&1"


# FIXME: single-end version needed
# NOTE: fastqc does not process reads in pairs. files are processed as single units.
# ANNOT: Do quality control on single end sample fastq read files.
rule fastqc_raw_se:
    input: trim_reads_input
    output:
        # all outputs are provided to ensure atomicity
        rep = os.path.join(FASTQC_DIR, '{sample}', '{sample}_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_fastqc.zip'),
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_raw.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    run:
        # renaming the ".fastq.gz" suffix to "_fastqc.html"
        tmp_output = os.path.basename(input[0]).replace(fastq_ext(input[0]), '_fastqc.html')
        tmp_zip = os.path.basename(input[0]).replace(fastq_ext(input[0]), '_fastqc.zip')
        shell("""{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1;
                if [[ {tmp_output} != {wildcards.sample}_fastqc.html ]]; then
                    mv {params.output_dir}/{tmp_output} {output.rep} &&\
                    mv {params.output_dir}/{tmp_zip} {output.zip}
                fi """)


# FIXME: or discard completely and change multiqc to use fastp --> fastp rule would have to be adjusted to create reasonable outputs
# ANNOT: Do quality control on paired end sample fastq read files.
# FIXME: Should this rule be called `fastqc_raw_pe`?
rule fastqc_raw:
    input: trim_reads_input
    output:
        r1_rep = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R1_fastqc.html'),
        r1_zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R1_fastqc.zip'),
        r2_rep = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R2_fastqc.html'),
        r2_zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R2_fastqc.zip') # all outputs are provided to ensure atomicity
    log: [os.path.join(LOG_DIR, 'fastqc_{sample}_raw_R1.log'), os.path.join(LOG_DIR, 'fastqc_{sample}_raw_R2.log')]
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    run:
        # renaming the ".fastq.gz" suffix to "_fastqc.html"
        tmp_R1_output = os.path.basename(input[0]).replace(fastq_ext(input[0]), '_fastqc.html')
        tmp_R1_zip = os.path.basename(input[0]).replace(fastq_ext(input[0]),  '_fastqc.zip')
        tmp_R2_output = os.path.basename(input[1]).replace(fastq_ext(input[0]),'_fastqc.html')
        tmp_R2_zip = os.path.basename(input[1]).replace(fastq_ext(input[0]),'_fastqc.zip')
        shell("""{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1;
                if [[ {tmp_R1_output} != {wildcards.sample}_R1_fastqc.html ]]; then
                    mv {params.output_dir}/{tmp_R1_output} {output.r1_rep} &&\
                    mv {params.output_dir}/{tmp_R1_zip} {output.r1_zip} &&\
                    mv {params.output_dir}/{tmp_R2_output} {output.r2_rep} &&\
                    mv {params.output_dir}/{tmp_R2_zip} {output.r2_zip}
                fi """)

# TODO: can probably be done by using map_input, no seperate functions neccessary?
rule fastqc_trimmed_se:
    input: os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed.fastq.gz")
    output:
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_fastqc.zip')
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_trimmed.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    shell: "{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1"

rule fastqc_trimmed_pe:
    input: os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R{read_num}.fastq.gz")
    output:
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_R{read_num}_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_R{read_num}_fastqc.zip')
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_trimmed_R{read_num}.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    shell: "{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1"

rule fastqc_primer_trimmed:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bam')
    output:
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.zip'),
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_aligned_primer-trimmed.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    shell: "{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1"

# TODO think about adding a global version to include all samples
rule multiqc:
  input: multiqc_input
  output: os.path.join(MULTIQC_DIR, '{sample}', 'multiqc_report.html')
  params:
    output_dir = os.path.join(MULTIQC_DIR, '{sample}')
  log: os.path.join(LOG_DIR, 'multiqc_{sample}.log')
  shell: "{MULTIQC_EXEC} -f -o {params.output_dir} {input} >> {log} 2>&1"


# TODO it should be possible to add customized parameter
rule lofreq:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bai'),
        ref = os.path.join(INDEX_DIR, "{}".format(os.path.basename(REFERENCE_FASTA)))
    output: vcf = os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    log: os.path.join(LOG_DIR, 'lofreq_{sample}.log')
    run:
        shell("{LOFREQ_EXEC} call -f {input.ref} -o {output} --verbose {input.aligned_bam} >> {log} 2>&1")
        # WIP create a dummy entry if no variant is found - use this as long as the input-function solution doesn't work
        no_variant_vep(wildcards.sample, output.vcf)


rule vcf2csv:
    input: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.csv')
    params:
        script = os.path.join(SCRIPTS_DIR, 'vcfTocsv.py')
    log: os.path.join(LOG_DIR, 'vcf2csv_{sample}.log')
    shell: "{PYTHON_EXEC} {params.script} {input} >> {log} 2>&1"

rule vep:
    input: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    output: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2.txt')
    params:
        species = "sars_cov_2"
    log: os.path.join(LOG_DIR, 'vep_{sample}.log')
    shell:
      """
      {VEP_EXEC} --verbose --offline --dir_cache {VEP_DB} --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing\
      --distance 5000 --mane --protein --species {params.species} --symbol --transcript_version --tsl\
      --input_file {input} --output_file {output} >> {log} 2>&1
      """

rule parse_vep:
    input: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2.txt')
    output: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2_parsed.txt')
    params:
        script = os.path.join(SCRIPTS_DIR, 'parse_vep.py')
    log: os.path.join(LOG_DIR, 'parse_vep_{sample}.log')
    shell: "{PYTHON_EXEC} {params.script} {input} {output} >> {log} 2>&1"


rule bam2fastq:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq')
    log: os.path.join(LOG_DIR, 'bam2fastq_{sample}.log')
    shell: "{SAMTOOLS_EXEC} fastq {input} > {output} 2>> {log} 3>&2"


rule kraken:
    input:
        unaligned_fastq = os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq'),
        database = KRAKEN_DB
    output: os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt')
    log: os.path.join(LOG_DIR, 'kraken_{sample}.log')
    shell: "{KRAKEN2_EXEC} --report {output} --db {input.database} {input.unaligned_fastq} >> {log} 2>&1"


rule krona_report:
    input:
        kraken_output = os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt'),
        database = KRONA_DB
    output: os.path.join(REPORT_DIR, '{sample}.Krona_report.html')
    log: os.path.join(LOG_DIR, 'krona_report_{sample}.log')
    shell: "{IMPORT_TAXONOMY_EXEC} -m 3 -t 5 {input.kraken_output} -tax {input.database} -o {output} >> {log} 2>&1"


rule samtools_bedcov:
    input:
        mutations_bed = MUTATIONS_BED,
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_mut_cov.csv')
    log: os.path.join(LOG_DIR, 'samtools_bedcov_{sample}.log')
    shell: "{SAMTOOLS_EXEC} bedcov {input.mutations_bed} {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule samtools_coverage:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_genome_cov.csv')
    log: os.path.join(LOG_DIR, 'samtools_coverage_{sample}.log')
    shell: "{SAMTOOLS_EXEC} coverage {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule create_sample_quality_table:
    input:
        genome_cov_file=os.path.join(COVERAGE_DIR, "{sample}_genome_cov.csv"),
        mut_cov_file=os.path.join(COVERAGE_DIR, "{sample}_mut_cov.csv"),
    output:
        os.path.join(COVERAGE_DIR, "{sample}_quality.csv"),
    params:
        script=os.path.join(SCRIPTS_DIR, "create_sample_quality_table.R"),
    log:
        os.path.join(LOG_DIR, "create_sample_quality_table_{sample}.log"),
    shell:
        """
        {RSCRIPT_EXEC} {params.script} \
            {wildcards.sample} \
            {input.genome_cov_file} \
            {input.mut_cov_file} \
            {output} > {log} 2>&1
        """


rule create_sample_quality_summary:
    input:
        script=os.path.join(SCRIPTS_DIR, "create_summary_table.R"),
        files=expand(
            os.path.join(COVERAGE_DIR, "{sample}_quality.csv"),
            sample=SAMPLES,
        ),
    output:
        os.path.join(COVERAGE_DIR, "sample_quality_summary.csv"),
    log:
        os.path.join(LOG_DIR, "create_sample_quality_summary.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} {output} {input.files} > {log} 2>&1
        """


rule generate_navbar:
    input:
      script = os.path.join(SCRIPTS_DIR, "generateNavigation.R")
    output:
      os.path.join(REPORT_DIR, "_navbar.html")
    params:
      report_scripts_dir = os.path.join(SCRIPTS_DIR, "report_scripts")
    log: os.path.join(LOG_DIR, "generate_navigation.log")
    shell: "{RSCRIPT_EXEC} {input.script} \
{params.report_scripts_dir} {SAMPLE_SHEET_CSV} {output} > {log} 2>&1"


rule render_kraken2_report:
    input:
      script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
      report=os.path.join(SCRIPTS_DIR, "report_scripts", "taxonomic_classification.Rmd"),
      header=os.path.join(REPORT_DIR, "_navbar.html"),
      kraken=os.path.join(KRAKEN_DIR, "{sample}_classified_unaligned_reads.txt"),
      krona=os.path.join(REPORT_DIR, "{sample}.Krona_report.html")
    output: os.path.join(REPORT_DIR, "{sample}.taxonomic_classification.html")
    log: os.path.join(LOG_DIR, "reports", "{sample}_taxonomic_classification.log")
    shell: """{RSCRIPT_EXEC} {input.script} \
{input.report} {output} {input.header} \
'{{\
  "sample_name": "{wildcards.sample}",  \
  "site_dir":    "{REPORT_DIR}",        \
  "krona_file":  "{input.krona}",       \
  "kraken_file": "{input.kraken}",      \
  "logo": "{LOGO}" \
}}' > {log} 2>&1"""


rule run_deconvolution:
    input:
        script=os.path.join(SCRIPTS_DIR, "deconvolution.R"),
        deconvolution_functions=os.path.join(
            SCRIPTS_DIR, "deconvolution_funs.R"
        ),
        vep=os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2_parsed.txt"),
        snv=os.path.join(VARIANTS_DIR, "{sample}_snv.csv"),
    output:
        sigmut_df=os.path.join(MUTATIONS_DIR, "{sample}_sigmuts.csv"),
        non_sigmut_df=os.path.join(MUTATIONS_DIR, "{sample}_non_sigmuts.csv"),
        variant_proportions=os.path.join(
            VARIANTS_DIR, "{sample}_variant_abundance.csv"
        ),
        variants_with_meta=os.path.join(VARIANTS_DIR,
         "{sample}_variants_with_meta.csv"),
        mutations=os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"),
    log:
        os.path.join(LOG_DIR, "{sample}_deconvolution.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
        "{input.deconvolution_functions}" \
        "{wildcards.sample}" \
        "{MUTATION_SHEET_CSV}" \
        "{SAMPLE_SHEET_CSV}" \
        "{input.vep}" \
        "{input.snv}" \
        "{MUTATION_DEPTH_THRESHOLD}" \
        "{output.sigmut_df}" \
        "{output.non_sigmut_df}" \
        "{output.variant_proportions}" \
        "{output.variants_with_meta}" \
        "{output.mutations}" \
        > {log} 2>&1
        """


rule render_variant_report:
    input:
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "variantreport_p_sample.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
        sigmut_file=os.path.join(MUTATIONS_DIR, "{sample}_sigmuts.csv"),
        non_sigmut_file=os.path.join(MUTATIONS_DIR, "{sample}_non_sigmuts.csv"),
        variant_abundance_file=os.path.join(
            VARIANTS_DIR, "{sample}_variant_abundance.csv"
        ),
        mutations=os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"),
        vep=os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2_parsed.txt"),
        snv=os.path.join(VARIANTS_DIR, "{sample}_snv.csv"),
    output:
        varreport=os.path.join(REPORT_DIR, "{sample}.variantreport_p_sample.html"),
    log:
        os.path.join(LOG_DIR, "reports", "{sample}_variant_report.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
        {input.report} {output.varreport} {input.header} \
        '{{ \
          "sample_name": "{wildcards.sample}", \
          "sigmut_file": "{input.sigmut_file}", \
          "non_sigmut_file": "{input.non_sigmut_file}", \
          "variant_abundance_file": "{input.variant_abundance_file}", \
          "snv_file": "{input.snv}", \
          "vep_file": "{input.vep}", \
          "logo": "{LOGO}" \
        }}' > {log} 2>&1
        """


rule render_qc_report:
    input:
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "qc_report_per_sample.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
        coverage=os.path.join(COVERAGE_DIR, "{sample}_quality.csv"),
        multiqc=os.path.join(MULTIQC_DIR, "{sample}", "multiqc_report.html"),
    output:
        html_report=os.path.join(REPORT_DIR, "{sample}.qc_report_per_sample.html"),
        table_outfile=os.path.join(
            COVERAGE_DIR, "{sample}_report_download_coverage.csv"
        ),
    params:
        multiqc_rel_path=lambda wildcards, input: input.multiqc[len(REPORT_DIR) + 1 :],
    log:
        os.path.join(LOG_DIR, "reports", "{sample}_qc_report.log"),
    shell:
        """{RSCRIPT_EXEC} {input.script} \
        {input.report} {output.html_report} {input.header} \
        '{{\
          "sample_name": "{wildcards.sample}",  \
          "coverage_file": "{input.coverage}",   \
          "multiqc_report": "{params.multiqc_rel_path}", \
          "logo": "{LOGO}", \
          "coverage_table_outfile": "{output.table_outfile}" \
        }}' > {log} 2>&1"""


rule create_variants_summary:
    input:
        script = os.path.join(SCRIPTS_DIR, "create_summary_table.R"),
        files = expand(os.path.join(VARIANTS_DIR, "{sample}_variants_with_meta.csv"), sample = SAMPLES)
    output: os.path.join(VARIANTS_DIR, 'data_variant_plot.csv')
    log: os.path.join(LOG_DIR, "create_variants_summary.log")
    shell: """
        {RSCRIPT_EXEC} {input.script} {output} {input.files} > {log} 2>&1
        """

rule create_mutations_summary:
    input:
        script = os.path.join(SCRIPTS_DIR, "create_summary_table.R"),
        files = expand(os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"), sample = SAMPLES)
    output: os.path.join(MUTATIONS_DIR, 'data_mutation_plot.csv')
    log: os.path.join(LOG_DIR, "create_mutations_summary.log")
    shell: """
        {RSCRIPT_EXEC} {input.script} {output} {input.files} > {log} 2>&1
        """


# TODO integrate the output of fastp.json to get the number of raw and trimmed reads
rule create_overviewQC_table:
    input:
        script=os.path.join(SCRIPTS_DIR, "overview_QC_table.R"),
        quality_table_file=os.path.join(COVERAGE_DIR, "sample_quality_summary.csv"),
    output:
        os.path.join(OUTPUT_DIR, "overview_QC.csv"),
    log:
        os.path.join(LOG_DIR, "create_overviewQC_table.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
            {SAMPLE_SHEET_CSV} \
            {output} \
            {READS_DIR} \
            {TRIMMED_READS_DIR} \
            {MAPPED_READS_DIR} \
            {input.quality_table_file} > {log} 2>&1
        """


rule run_mutation_regression:
    input:
        script=os.path.join(SCRIPTS_DIR, "mutation_regression.R"),
        mutations_csv=os.path.join(MUTATIONS_DIR, "data_mutation_plot.csv"),
        overviewQC=os.path.join(OUTPUT_DIR, "overview_QC.csv"),
        fun_lm=os.path.join(SCRIPTS_DIR, "pred_mutation_increase.R"),
        fun_tbls=os.path.join(SCRIPTS_DIR, "table_extraction.R"),
        mutation_sheet=MUTATION_SHEET_CSV
    output:
        mut_count_outfile=os.path.join(OUTPUT_DIR, "mutations_counts.csv"),
        unfilt_mutation_sig_outfile=os.path.join(
            OUTPUT_DIR, "unfiltered_mutations_sig.csv"
        ),
    log:
        os.path.join(LOG_DIR, "mutation_regression.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
            {input.mutations_csv} \
            {input.mutation_sheet} \
            {input.fun_lm} \
            {input.fun_tbls} \
            {MUTATION_COVERAGE_THRESHOLD} \
            {input.overviewQC} \
            {output.mut_count_outfile} \
            {output.unfilt_mutation_sig_outfile} \
            > {log} 2>&1
        """


rule render_index:
    input:
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "index.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
        variants=os.path.join(VARIANTS_DIR, "data_variant_plot.csv"),
        mutations=os.path.join(MUTATIONS_DIR, "data_mutation_plot.csv"),
        overviewQC=os.path.join(OUTPUT_DIR, "overview_QC.csv"),
        mut_count_file=os.path.join(OUTPUT_DIR, "mutations_counts.csv"),
        unfiltered_mutation_sig_file=os.path.join(
            OUTPUT_DIR, "unfiltered_mutations_sig.csv"
        ),
    params:
        fun_tbls=os.path.join(SCRIPTS_DIR, "table_extraction.R"),
        fun_pool=os.path.join(SCRIPTS_DIR, "pooling.R"),
        fun_index=os.path.join(SCRIPTS_DIR, "fun_index.R"),
    output:
        report=os.path.join(REPORT_DIR, "index.html"),
    log:
        os.path.join(LOG_DIR, "reports", "index.log"),
    shell:
        """{RSCRIPT_EXEC} {input.script} \
        {input.report} {output.report} {input.header}   \
        '{{ \
          "variants_csv": "{input.variants}", \
          "mutations_csv": "{input.mutations}", \
          "sample_sheet": "{SAMPLE_SHEET_CSV}", \
          "mutation_sheet": "{MUTATION_SHEET_CSV}", \
          "mutation_coverage_threshold": "{MUTATION_COVERAGE_THRESHOLD}", \
          "mut_count_file": "{input.mut_count_file}", \
          "unfiltered_mutation_sig_file": "{input.unfiltered_mutation_sig_file}", \
          "logo": "{LOGO}", \
          "fun_tbls": "{params.fun_tbls}", \
          "fun_pool": "{params.fun_pool}", \
          "fun_index": "{params.fun_index}", \
          "overviewQC": "{input.overviewQC}", \
          "output_dir": "{OUTPUT_DIR}" \
        }}' > {log} 2>&1"""
