import subprocess
import numpy
import os

def figshare_download(save_path):
    process = subprocess.run(['wget', '-O', save_path,
                    'https://figshare.com/ndownloader/files/44999293'], check=True)
    # If the downloaded filename ends in tar.gz then extraact it
    if save_path.endswith(".tar.gz"):
       with tarfile.open(save_path) as tar:
            tar.extractall(path=os.path.dirname(save_path))
            print("Done!")