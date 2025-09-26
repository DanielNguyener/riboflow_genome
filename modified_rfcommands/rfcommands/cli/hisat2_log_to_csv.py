# -*- coding: utf-8 -*-

from .main import *
from ..hisat2_log_to_csv import parse_hisat2_log, write_csv
import sys

@cli.command()
@click.option('--log', '-l',
              type     = click.Path(exists = True),
              help     = "HISAT2 log file.",
              required = True)
@click.option('--name', '-n',
              type     = click.STRING,
              help     = "Name of the experiment (column label).",
              required = True )
@click.option('--prefix', '-p',
              type     = click.STRING,
              help     = "Prefix for row names.",
              required = True )
@click.option('--out', '-o',
              help     = "Output CSV file.",
              type     = click.Path(exists = False),
              required = True)
def hisat2_log_to_csv(log, name, prefix, out):
    """
    Converts HISAT2 alignment statistics to CSV format.
    """
    try:
        stats = parse_hisat2_log(log)
        write_csv(stats, out, name, prefix)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)