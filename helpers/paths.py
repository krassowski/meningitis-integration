from pathlib import Path


def make_directories(dirs, root):
    root_path = Path(root)
    for directory, subdirectories in dirs.items():
        directory_path = root_path / directory
        for subdirectory in subdirectories:
            path = directory_path / subdirectory
            path.mkdir(exist_ok=True)
