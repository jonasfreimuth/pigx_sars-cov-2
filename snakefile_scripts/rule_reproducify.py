import json
import pathlib
import sys
import os
import shutil
import subprocess
from snakemake.logging import logger

with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = sys.stdout    

    target_files = dict(zip(
        snakemake.params["output_keys"],
        snakemake.output))

    config = snakemake.params["config"]

    # Write the config file
    with open(target_files["config_out"], 'w') as outfile:
        dumps = json.dumps(config,
                           indent=4, sort_keys=True,
                           separators=(",",": "), ensure_ascii=True)
        outfile.write(dumps)

    # Copy "parameter files" into reproducify dir
    param_file_keys = snakemake.params["param_file_keys"]

    for key in param_file_keys:
        origin = config["locations"][key]
        target = target_files[key]
        shutil.copy2(origin, target)

    # get the pipeline version string
    # This is very roundabout owing to the fact that neither the location of the
    # pipeline executable, nor the VERSION file have fixed locations.
    src_dir = config["locations"]["exec_prefix"]
    src_sub = "bin/"

    if os.getenv("PIGX_UNINSTALLED"):
        src_sub = ""

    pipeline_version_call = [
        f"{src_dir}/{src_sub}pigx-sars-cov-2",
        f"--version"]
    
    logger.info(pipeline_version_call)

    try:
        with open(target_files["pigx_version_out"], "w") as version_out_file:
            subprocess.run(
                pipeline_version_call,
                stdout = version_out_file,
                stderr = log_file)
    except Exception as e:
        logger.info(e)
        sys.exit(1)

    # TODO Add database versions
    for file in snakemake.output:
        filepath = pathlib.Path(file)
        filepath.touch()