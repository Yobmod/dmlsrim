
from git import Repo
from pathlib import Path
import dotenv
import os
from datetime import datetime
import subprocess
import json

# TODO: get username and pw using dotenv
# use typer or argparse to get commit msg. if blank, use date + id?
# Do a check with "safety" lib. use subprocess?
# add progressbar from git lib?

cwd = Path.cwd().resolve()
repo = Repo(cwd)
assert not repo.bare
assert repo.working_tree_dir == str(cwd.resolve())

dotenv.load_dotenv(cwd / ".env")
username = os.getenv("GIT_USERNAME") or "DML"
password = os.getenv("GIT_PASSWORD")  # or None

try:
    repo.remotes.origin.pull()
    print("Pulling from github....")
except Exception as e:
    print("Could not pull from github!")
    raise e
else:
    print("                  .... done")

if repo.is_dirty():
    changedFiles = [item.a_path for item in repo.index.diff(None)]
    repo.git.add(update=True)
    print(f"Modifying files: {changedFiles}")
    # TODO get list of changed files. repo.diff? repo.git.diff?
    if repo.untracked_files:
        print(f"Adding files: {repo.untracked_files}")  # list of filename strings that have not been added
        repo.index.add(repo.untracked_files)
    update_pending = True
else:
    update_pending = False


commit_msg = ""


if update_pending:
    if not commit_msg:
        commit_msg = f"{username}, {datetime.now()}"
    try:
        repo.index.commit(commit_msg)
        print(f"Commiting changes with msg: {commit_msg}")
        repo.remotes.origin.push()
        print(f"Pushing changes to: {repo.remotes.origin.url}")
    except Exception:
        raise

venv_bytespath = subprocess.check_output("poetry env info --path".split(), shell=True)
venv_str = venv_bytespath.decode("UTF-8")
venv_path = venv_str.strip()

# check if vscode settings exists, or create if not
Path(".vscode").mkdir(parents=True, exist_ok=True)
Path(".vscode/settings.json").touch()

print(f"Updating path to \n {venv_path}")

with open(".vscode/settings.json", "r") as f:
    settings = json.load(f)
    settings["python.pythonPath"] = venv_path + R"\\scripts\\python.exe"

with open(".vscode/settings.json", "w") as f:
    json.dump(settings, f, sort_keys=True, indent=4)


# print(json.dumps(settings, sort_keys=True, indent=4))
