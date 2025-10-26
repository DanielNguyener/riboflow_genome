# -*- coding: utf-8 -*-

from .main import *

from ..compile_step_stats import compile

@cli.command()
@click.option('--out', '-o',
              help     = "Output file. The stats are written in csv format"
                         "to this file." ,
              required = True,
              type     = click.Path(exists = False))
@click.option('--cutadapt', '-c',
              help     = "Cutadapt log file." ,
              required = True,
              type     = click.Path(exists = True))
@click.option('--filter', '-f',
              help     = "Filter alignment log file." ,
              required = True,
              type     = click.Path(exists = True))          
@click.option('--trans', '-t',
              help     = "Transcriptome alignment log file." ,
              required = False,
              type     = click.Path(exists = True))
@click.option('--genome', '-g',
              help     = "Genome alignment log file." ,
              required = False,
              type     = click.Path(exists = True))
@click.option('--quality', '-q',
              help     = "Quality/passing reads count file.",
              required = True,
              type     = click.Path(exists = True))
@click.option('--dedup', '-d',
              help     = "Deduplication count file. "
                         "This file should only contain one number." ,
              required = True,
              type     = click.Path(exists = True))
@click.option('--psite', '-p',
              help     = "P-site count file. This file should only contain one number." ,
              required = False,
              type     = click.Path(exists = True))
@click.option('--label-prefix', '-l',
              help     = "Label prefix for statistics (transcriptome or genome)." ,
              default = "transcriptome",
              type     = click.STRING)
@click.option('--name', '-n',
              help     = "Name of the experiment." ,
              required = True,
              type     = click.STRING)
def compile_step_stats(out,     cutadapt, filter, trans, genome,
                       quality,   dedup,  psite, label_prefix, name):
    """
    Puts statistics coming from various steps into one file.

    Merges cutadapt, and alignment statistics coming from Bowtie2 or Hisat
    from the given samples into one file.
    This is done by summing up the corresponding counts
    and calculating the percentages.
    This version is implemented for single-end reads only.
    So it won't work for paired-end statistics yet.

    For convenience, we are providing percentages of these statistics.
    We are rounding up all percentages to integers for simplicity.
    If you want higher precision,
    you can re-calculate using the counts given in these tables.
    """

    compile(out      = out,
            cutadapt = cutadapt,
            filter   = filter,
            trans    = trans,
            genome   = genome,
            quality  = quality,
            dedup    = dedup,
            psite    = psite,
            label_prefix = label_prefix,
            name     = name)
