import typing as t
from typer import Typer, echo

app = Typer()


@app.command()
def hello(name: str) -> None:
    echo(f"Hello {name}")


@app.command()
def goodbye(name: str, formal: bool = False) -> None:
    if formal:
        echo(f"Goodbye Ms. {name}. Have a good day.")
    else:
        echo(f"Bye {name}!")


@app.command()
def run_dcf(filename: str, formal: bool = False) -> None:
    if formal:
        echo(f"Goodbye Ms. {filename}. Have a good day.")
    else:
        echo(f"Bye {filename}!")


if __name__ == "__main__":
    app()
