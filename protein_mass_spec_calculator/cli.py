"""Custom CLI for protein-mass-spec-calculator.

This is totally optional;
if you want to use it, though,
follow the skeleton to flesh out the CLI to your liking!
Finally, familiarize yourself with Typer,
which is the package that we use to enable this magic.
Typer's docs can be found at:

    https://typer.tiangolo.com
"""

import typer

app = typer.Typer()


@app.command()
def hello():
    """Echo the project's name."""
    typer.echo("This project's name is protein-mass-spec-calculator")


@app.command()
def describe():
    """Describe the project."""
    typer.echo("Omni calculator for proteins + mass spec")


if __name__ == "__main__":
    app()
