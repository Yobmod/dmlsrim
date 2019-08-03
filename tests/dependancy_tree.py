
import contextlib
import sys
from io import StringIO
from pathlib import Path

from typing import ContextManager, Generator


@contextlib.contextmanager
def captured_output(stream_name: str) -> Generator[StringIO, None, None]:
    """Return a context manager used by captured_stdout and captured_stdin
    that temporarily replaces the sys stream *stream_name* with a StringIO."""
    orig_stdout = getattr(sys, stream_name)
    setattr(sys, stream_name, StringIO())
    try:
        yield getattr(sys, stream_name)
    finally:
        setattr(sys, stream_name, orig_stdout)


def captured_stdout() -> ContextManager[StringIO]:
    """Capture the output of sys.stdout"""
    return captured_output("stdout")


def captured_stderr() -> ContextManager[StringIO]:
    """Capture the output of sys.strerr"""
    return captured_output("stderr")


def captured_stdin() -> ContextManager[StringIO]:
    """Capture the output of sys.stdin"""
    return captured_output("stdin")


# have to get cwd before import pipenv
cwd = Path(__file__).absolute()
save_path = cwd.parent / 'graph.txt'
json_path = cwd.parent / 'graph.json'


try:
    save_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)
except Exception as e:
    print(e)


with open(save_path, 'w') as file:
    with captured_stdout() as stream:
        try:
            from pipenv.core import do_graph    # NOQA
            do_graph(bare=False, json=False, json_tree=False, reverse=True)
        except SystemExit:   # core.do_graph() calls sys.exit()
            str_out = stream.getvalue()
            file.write(str_out)
    print(f'Dependancy graph captured to {save_path}')


"""
with open(json_path, 'w') as file:
    with captured_stdout() as stream:
        try:
            from pipenv.core import do_graph    # NOQUA
            do_graph(bare=False, json=False, json_tree=True, reverse=False)
        except SystemExit:   # core.do_graph() calls sys.exit()
            str_out = stream.getvalue()
            file.write(str_out)
    print(f'Dependancy graph captured to {json_path}')
"""
