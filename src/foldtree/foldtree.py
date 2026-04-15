#!/usr/bin/env python3

"""
Foldtree creates phylogenetic trees from protein 
structures using Foldseek to align protein structures 
and generate distance matrices for tree construction.

This file wraps the snakemake pipeline to make it
easier to call from the command line.
"""

import click
import os
import shlex
import subprocess
from importlib.metadata import version as pkg_version, PackageNotFoundError
from pathlib import Path
from typing import Optional


def get_version() -> str:
    try:
        return pkg_version("foldtree")
    except PackageNotFoundError:
        return "0+unknown"    


def _bool_to_str(v: bool) -> str:
    return "True" if v else "False"


def _get_snakemake() -> Path:
    """Returns snakemake path"""
    conda_prefix = os.environ.get("CONDA_PREFIX")
    p = Path(conda_prefix) / "bin" / "snakemake"
    return p

def _choose_snakefile() -> Path:
    snakefile = _snakefile_from_sources()
    return _validate_snakefile(snakefile)


def _snakefile_from_sources() -> Path:
    """
    Returns the path of the main Foldtree snakefile
    from this directory; that is, the one shipped with this
    exact file. To run foldtree from sources
    """
    snakefile = Path(f"{Path(__file__).resolve().parent}/workflow/fold_tree")
    return Path(snakefile)



def _validate_snakefile(p: Path) -> Path:
    if not p.exists():
        raise click.ClickException(f"Snakefile not found: {p}")
    if p.is_dir():
        raise click.ClickException(f"Snakefile path is a directory, expected a file: {p}")
    return p


def _check_snakemake() -> None:
    try:
        snakemake = _get_snakemake()
        subprocess.run([snakemake, "-h"], 
                       check=True, 
                       stderr=subprocess.DEVNULL, 
                       stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        raise click.ClickException(
            "Could not find 'snakemake' on PATH. Something went wrong during installation."
        )

def _validate_folder(folder: Path) -> Path:
    if not folder.exists():
        raise click.ClickException(f"Folder does not exist: {folder}")
    if not folder.is_dir():
        raise click.ClickException(f"Folder is not a directory: {folder}")
    return folder


# Foldtree-related options
@click.command(context_settings={"help_option_names": ["-h", "--help"]},
               no_args_is_help=True)
@click.option(
    "--folder",
    "folder",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True, exists=True),
    required=True,
    help="Input folder with identifiers.txt and optional 'structs' folder.",
)
@click.option(
    "--filter/--no-filter",
    default=False,
    show_default=True,
    help="Filter out structures with an average PLDDT < 40 (improves tree quality, see manuscript).",
)
@click.option(
    "--foldseek-cores",
    type=int,
    default=None,
    help="Number of parallel jobs for all-vs-all comparison by Foldseek.",
)
@click.option(
    "--custom-structs", "-cs",
    is_flag=True,
    default=False,
    help="Use a precompiled sets of structures defined in FOLDER/structs. FOLDER/identifiers.txt should be empty.",
)
# Snakemake-related options
@click.option(
    "--workdir",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True),
    default=None,
    help="Run snakemake from this working directory (defaults to current directory).",
)
@click.option(
    "--cores", "-c",
    type=int,
    default=4,
    show_default=True,
    help="Number of parallel jobs run by snakemake.",
)
@click.option(
    "--dry-run", "-n",
    is_flag=True,
    default=False,
    help="Run snakemake with -n/--dry-run. Does not execute anything, and displays what would be done.",
)
@click.option(
    "--printshellcmds", "-p",
    is_flag=True,
    default=False,
    help="Run snakemake with -p/--printshellcmds. Prints commands run by snakemake.",
)
@click.option(
    "--extra-snakemake-args", "-esa",
    default="",
    help="Extra args passed verbatim to snakemake (e.g. '--rerun-incomplete' or cluster profile settings) for custom configuration.",
)
# Other
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Print the snakemake command before running.",
)
@click.version_option(version=get_version())
def cli(
    folder: Path,
    filter: bool,
    foldseek_cores: Optional[int],
    custom_structs: bool,
    workdir: Optional[Path],
    cores: int,
    dry_run: bool,
    printshellcmds: bool,
    extra_snakemake_args: str,
    verbose: bool
) -> None:
    """
    Foldtree creates phylogenetic trees from protein
    structures using Foldseek to align protein structures
    and generate distance matrices for tree construction.

    \b
    Usage:
        foldtree --folder <input_folder> --cores <N>

    \b
    Example:
        # Default local run
        foldtree --folder myfam -p

    \b
        # Single thread
        foldtree --folder myfam -p -c 1

    \b
        # Run on a SLURM cluster:
        foldtree --folder myfam -p -c 8 -esa "--profile slurmsimple/simple"

    \b
    Please cite:
        Moi D, Bernard C, Steinegger M, Nevers Y, Langleib M, Dessimoz C.
        Structural phylogenetics unravels the evolutionary diversification
        of communication systems in gram-positive bacteria and their viruses.
        Nature Structural & Molecular Biology
    \b
        doi: 10.1038/s41594-025-01649-8 (2025)

    \b
    Contributors:
        David Moi, Mauricio Langleib, Martin Steinegger, 
        Stefano Pascarelli, Nikolai Romashchenko
        (c) 2022-present
    """
    if cores < 1:
        raise click.ClickException("--cores must be >= 1")

    if foldseek_cores is not None and foldseek_cores < 1:
        raise click.ClickException("--foldseek-cores must be >= 1")

    folder = _validate_folder(folder)

    snakefile = _choose_snakefile()
    print("Running pipeline:", snakefile)

    cmd = ["snakemake"]
    cmd += ["--cores", str(cores)]
    cmd += ["-s", str(snakefile)]

    # snakemake CLI config values are strings, so be explicit.
    config_kv = [
        f"folder={str(folder)}",
        f"filter={_bool_to_str(filter)}",
    ]
    if foldseek_cores is not None:
        config_kv.append(f"foldseek_cores={foldseek_cores}")
        
    if custom_structs:
        config_kv.append("custom_structs=True")

    cmd += ["--config", *config_kv]

    if dry_run:
        cmd.append("--dry-run")
    
    if printshellcmds:
        cmd.append("--printshellcmds")

    if extra_snakemake_args.strip():
        # Split the arguments similar to how shell would do it
        cmd += shlex.split(extra_snakemake_args)

    run_cwd = str(workdir) if workdir else None

    if verbose:
        click.echo("Running:\n  " + " ".join(shlex.quote(x) for x in cmd))
        if run_cwd:
            click.echo(f"Working directory: {run_cwd}")

    # check if snakemake is available
    _check_snakemake()

    print(" ".join(x for x in cmd))
    try:
        subprocess.run(cmd, check=True, cwd=run_cwd)
    except subprocess.CalledProcessError as e:
        raise click.ClickException(f"Snakemake failed with exit code {e.returncode}")


if __name__ == "__main__":
    cli()
