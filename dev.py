from git import Repo
from pathlib import Path
import dotenv
import os
from datetime import datetime
# import safety

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
