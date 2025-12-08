# -*- coding: utf-8 -*-

from .main import *
from ..merge_hisat2_logs import merge_stats, write_merged_log
import sys

@cli.command()
@click.option('--output', '-o',
              help     = "Output merged log file.",
              type     = click.Path(exists = False),
              required = True)
@click.argument('logs', nargs=-1, required=True)
def merge_hisat2_logs(output, logs):
    """
    Merge multiple HISAT2 log files into a single log file.
    """
    try:
        merged_stats = merge_stats(logs)
        write_merged_log(merged_stats, output)
        print(f"Merged {len(logs)} log files into {output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

