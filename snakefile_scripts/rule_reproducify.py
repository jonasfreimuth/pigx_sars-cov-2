from fileinput import filename
import json
import pathlib
import sys
import subprocess
from snakemake.logging import logger

with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = sys.stdout    

    target_files = dict(zip(
        snakemake.params["output_keys"],
        snakemake.output))
    
    logger.info(target_files)

    # TODO Actually fill these files
    for file in snakemake.output:
        filepath = pathlib.Path(file)
        filepath.touch()