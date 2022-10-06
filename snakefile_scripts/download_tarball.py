import re
import tarfile
import sys
import urllib.request as urlrequest
from snakemake.logging import logger


def download_tarball(dl_url):
    """
    Download .tar.gz archive from url and return a tarfile object to it.
    """

    if not re.search(".tar.gz$", dl_url):
        raise Exception(f"dl_url ({dl_url}) does not point to a '.tar.gz' file.")

    if not re.match("[a-z]*://", dl_url):
        logger.info(
            f"No protocol identifier in {dl_url}, assuming http/https...")
        dl_url = "http://" + dl_url

    logger.info(
        f"Attempting to download database archive from {dl_url}...")
    
    try:
        dl_resp = urlrequest.urlopen(dl_url)

    except Exception as e:
        logger.error(e)
        sys.exit(1)

    tar_out = tarfile.open(fileobj = dl_resp, mode="r:gz")

    return(tar_out)
