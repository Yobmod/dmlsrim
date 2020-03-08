from git import Repo
from pathlib import Path
import dotenv

# TODO: get username and pw using dotenv
# use typer or argparse to get commit msg. if blank, use date + id?


cwd = Path.cwd()

repo = Repo(cwd)
assert not repo.bare
#
# assert repo.working_tree_dir == cwd.resolve()

print(cwd.resolve())
print(repo.working_tree_dir)


"""provides high-level access to your data, 
it allows you to create and delete heads, 
tags and remotes and access the configuration of the repository"""

print(repo.untracked_files)  # list of filename strings that have not been added

repo.index.add(repo.untracked_files)
repo.index.commit("commit message")
repo.pull()
repo.push()
