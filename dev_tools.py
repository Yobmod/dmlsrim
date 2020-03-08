from git import Repo
from pathlib import Path
import dotenv
import os

# TODO: get username and pw using dotenv
# use typer or argparse to get commit msg. if blank, use date + id?


cwd = Path.cwd().resolve()
repo = Repo(cwd)
assert not repo.bare
assert repo.working_tree_dir == str(cwd.resolve())


dotenv.load_dotenv(cwd / ".env")
username = os.getenv("GIT_USERNAME")  # or None
password = os.getenv("GIT_PASSWORD")  # or None


if repo.is_dirty():
    repo.git.add(update=True)
    if repo.untracked_files:
        print(f"Adding files: {repo.untracked_files}")  # list of filename strings that have not been added
    update_pending = True
    # repo.index.add(repo.untracked_files)
else:
    update_pending = False

if update_pending:
    try:
        repo.index.commit("commit message2")
        repo.remotes.origin.push()
        print(f"Pushing changes to {repo.remotes.origin.url}")
    except Exception:
        raise
