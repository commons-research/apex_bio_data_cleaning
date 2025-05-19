import os
import subprocess
from typing import Literal

import environments_utils

from src.constants import LINUX_BIN_PATH, MACOS_BIN_PATH
from src.neg_batch_mode import generate_negative_batch_mode
from src.pos_batch_mode import generate_positive_batch_mode


def write_mzbatch_file(
    parent_dir: str,
    mzml_input_file: str,
    mzbatch_file_name: str,
    cas_number_ion_str: str,
    ionisation: Literal["pos", "neg"],
) -> None:
    if ionisation == "pos":
        batch_file = generate_positive_batch_mode(
            os.path.join(os.getcwd(), parent_dir, mzml_input_file),
            os.path.join(os.getcwd(), parent_dir, "feature_list.csv"),
            os.path.join(os.getcwd(), parent_dir, f"{cas_number_ion_str}_sirius.mgf"),
            os.path.join(os.getcwd(), parent_dir, f"{cas_number_ion_str}_gnps.mgf"),
            os.path.join(os.getcwd(), parent_dir, f"{cas_number_ion_str}.csv"),
        )
    elif ionisation == "neg":
        batch_file = generate_negative_batch_mode(
            os.path.join(os.getcwd(), parent_dir, mzml_input_file),
            os.path.join(os.getcwd(), parent_dir, "feature_list.csv"),
            os.path.join(os.getcwd(), parent_dir, f"{cas_number_ion_str}_sirius.mgf"),
            os.path.join(os.getcwd(), parent_dir, f"{cas_number_ion_str}_gnps.mgf"),
            os.path.join(os.getcwd(), parent_dir, f"{cas_number_ion_str}.csv"),
        )
    else:
        raise ValueError("Ionisation must be 'pos' or 'neg'.")

    with open(os.path.join(parent_dir, mzbatch_file_name), "w") as f:
        f.write(batch_file)


def run_mzmine(working_dir: str, mz_bacth_file: str):
    if environments_utils.is_macos():
        mzmine_bin_path = MACOS_BIN_PATH
    elif environments_utils.is_linux():
        mzmine_bin_path = LINUX_BIN_PATH
    else:
        raise ValueError(
            "For Windows, check the mzmine documentation for batch processing."
        )

    subprocess.run(
        [
            mzmine_bin_path,
            "-b",
            mz_bacth_file,
        ],
        cwd=working_dir,
    )
