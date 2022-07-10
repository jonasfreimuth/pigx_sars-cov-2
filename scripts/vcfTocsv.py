# 21/04/04
# Vic-Fabienne Schumann
# input: vcf file


import os
import re
import csv
import sys


def vcfTocsv(vcffile, csvfile):
    """
    This functions takes a vcf as input, parses the information
    given in it's "info" column and writes it to a given csv file so it can be
    loaded into an R file easily.
    It requires that the fields "DP", "AF", "SB", and "DP4" are present, which
    may not be the case for input files generated by programs other than lofreq.
    """

    # because i get the "has to have a write option error"
    with open(csvfile, "w", newline='') as csvfile:

        # DP4 is missing but don't know right now how to deal with it
        fieldnames = ["Chromosome", "Pos", "Ref",
                      "Var", "DP", "AF", "SB", "DP4"]

        writer = csv.DictWriter(csvfile,  # TODO: make this dynamic
                                fieldnames=fieldnames,
                                delimiter=",")
        writer.writeheader()

        with open(vcffile, "r") as vcf:

            for line in vcf:
                if re.match("^NC_", line):
                    snv_info1 = re.split(r'\t+', line)
                    Chrom = snv_info1[0]
                    Pos = snv_info1[1]
                    Ref = snv_info1[3]
                    Var = snv_info1[4]

                    snv_info2 = re.split(r';', snv_info1[7])
                    for element in snv_info2:
                        if re.search(r'DP=', element):
                            # can this be concat. to fewer lines?
                            DP = element.split("=")[1]
                        if re.search(r'AF=', element):
                            # this part feels to be a bit redundant
                            AF = element.split("=")[1]
                        if re.search(r'SB=', element):
                            SB = element.split("=")[1]
                        if re.search(r'DP4=', element):
                            DP4 = element.split("=")[1]

                    writer.writerow({'Chromosome': Chrom,
                                     'Pos': Pos,
                                     'Ref': Ref,
                                     'Var': Var,
                                     'DP': DP,
                                     'AF': AF,
                                     'SB': SB,
                                     # has exclamations in the output, maybe get
                                     # rid of them, don't know how right now,
                                     # also not super important right now
                                     'DP4': DP4.split("\n")[0]})


if __name__ == '__main__':
    vcffile = sys.argv[1]
    csvfile = sys.argv[2]
    vcfTocsv(vcffile, csvfile)
