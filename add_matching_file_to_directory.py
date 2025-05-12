
import shutil
import os 
path = os.getcwd()
directories = os.listdir(os.path.join(path, "example_folder"))

all_files = os.listdir(
    os.path.join(
        path,
        "data_files"
    )
)

for dir in directories:
    files = [file for file in all_files if dir in file] 

    if files:
        m = max(files, key=len) 
        file_path = os.path.join(path, "data_files", m))
        dest_dir = os.path.join(path, "example_folder", dir)
        shutil.move(file_path, dest_dir)

