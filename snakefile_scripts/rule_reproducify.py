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

    target_files = snakemake.params["targets"] 

    config = snakemake.params["config"]

    # Write the config file
    with open(target_files["config_out"], 'w') as outfile:
        dumps = json.dumps(config,
                           indent=4, sort_keys=True,
                           separators=(",",": "), ensure_ascii=True)
        outfile.write(dumps)

    # Copy "parameter files" into reproducify dir
    param_file_keys = snakemake.params["param_file_keys"]

    # Handle the special case of the db_dl.log file, which may not be there.
    if config["parameters"]["reproducify"]["save-db-version-log"]:
        db_dl_log_file = config["locations"]["db-dl-log-file"]

        if not os.path.exists(db_dl_log_file):
            script_file = (
                config["locations"]["pkglibexecdir"] +
                "/scripts/download_databases.sh")
            logger.error(
                f"ERROR: Pipeline is configured to save the log file of "
                f"database downloads to the\nreproducify dir, but could not "
                f"find the log file at\n{db_dl_log_file}.\nThe "
                f"file should be created automatically upon\n"
                f"dowloading the databases using the script at\n"
                f"{script_file}.\n"
                f"To disable this message set\n\"parameters: reproducify: "
                f"save-db-version-log\" to \"no\" in the settings file.")
            sys.exit(1)

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
        logger.error(e)
        sys.exit(1)

    for file in snakemake.output:
        filepath = pathlib.Path(file)
        filepath.touch()