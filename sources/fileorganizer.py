#!/usr/bin/env python3

import sys
from pathlib import Path


def organize_gz_files(base_path):
    if not base_path.is_dir():
        print(f"ERROR: The specified folder does not exist:\n{base_path}")
        sys.exit(1)

    print(f"Scanning: {base_path}\n")

    files_to_move = list(base_path.glob('*.gz'))

    if not files_to_move:
        print("No .gz files found to organize. Exiting.")
        return

    print(f"Found {len(files_to_move)} .gz files. Starting organization...")

    for file_path in files_to_move:
        new_folder_name = file_path.name.split('.')[0]

        if not new_folder_name:
            print(f"  [!] Skipping file with invalid name: {file_path.name}")
            continue

        destination_folder = base_path / new_folder_name
        new_file_path = destination_folder / file_path.name

        destination_folder.mkdir(exist_ok=True)

        try:
            file_path.rename(new_file_path)
            print(f"  [OK] Moved: {file_path.name:30s} ->  {new_folder_name}/")
        except Exception as e:
            print(f"  [!] ERROR while moving {file_path.name}: {e}")

    print("\nOrganization completed.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Error: You must provide the folder path as an argument.")
        print(f"Usage: python {sys.argv[0]} /path/to/your/folder")
        sys.exit(1)

    input_path = Path(sys.argv[1])
    organize_gz_files(input_path)
