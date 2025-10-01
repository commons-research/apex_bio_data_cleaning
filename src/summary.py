import os
from pathlib import Path

import pandas as pd


def is_quant_file_empty(filepath):
    """Return True if the file is empty (no rows), False otherwise."""
    try:
        df = pd.read_csv(filepath)
        return df.empty
    except pd.errors.EmptyDataError:
        return True
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return True


def check_column_exists(filepath, column_name):
    """Check if a specific column exists in the CSV file."""
    try:
        df = pd.read_csv(filepath, nrows=1)  # Only need headers
        return column_name in df.columns
    except Exception as e:
        print(f"Error checking columns in {filepath}: {e}")
        return False


def write_summary(summary, filename):
    if summary:
        df = pd.DataFrame(summary)
        df.to_csv(filename, index=False)
        print(f"Wrote summary to {filename}")
    else:
        print(f"No data to write for {filename}")
