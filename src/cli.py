import typer
import typing as t

app = typer.Typer()


@app.command()
def hello(name: str) -> None:
    typer.echo(f"Hello {name}")


@app.command()
def goodbye(name: str, formal: bool = False) -> None:
    if formal:
        typer.echo(f"Goodbye Ms. {name}. Have a good day.")
    else:
        typer.echo(f"Bye {name}!")


@app.command()
def run_dcf(filename: str, formal: bool = False) -> None:
    if formal:
        typer.echo(f"Goodbye Ms. {filename}. Have a good day.")
    else:
        typer.echo(f"Bye {filename}!")


if __name__ == "__main__":
    app()
