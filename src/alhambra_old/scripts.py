import argparse
import logging
import os.path
import sys

import click
import ruamel.yaml as yaml

from .tilesets import TileSet

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

POOL = None


@click.group(context_settings=CONTEXT_SETTINGS, chain=True)
@click.argument("filename", type=click.Path(exists=True))
@click.option("--verbose/--quiet", "-v/-q", default=False)
@click.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    default=None,
    help="output file name (defaults to stdout)",
)
def alhambra(filename, verbose, output):
    pass


@alhambra.resultcallback()
def process_pipeline(processors, filename, verbose, output):
    if verbose:
        logging.basicConfig(stream=sys.stderr, level=logging.INFO)
    ts = TileSet.from_file(filename)
    for processor in processors:
        ts = processor(ts)
    if output is None:
        ts.to_file(sys.stdout)
    else:
        ts.to_file(output)


@alhambra.command(short_help="generate sticky end sequences")
@click.option(
    "--method",
    "-m",
    type=click.Choice(["default", "multimodel"]),
    help="stickydesign method to use (multimodel or default)",
    default="default",
)
@click.option("--trials", "-n", default=100, help="number of trials (defaults to 100)")
def create_end_sequences(method, trials):
    """generates sticky end sequences for the system."""

    def processor(ts):
        ts_new, _ = ts.create_end_sequences(method=method, trials=trials)
        return ts_new

    return processor


@alhambra.command(short_help="attempt to reduce the # of sticky ends")
@click.option("--latticedefects/--no-latticedefects", "-l/-L", default=True)
@click.option("--sclass", "-s", multiple=True, default=["2GO"])
@click.option("--trials", "-n", default=1)
def reduce_ends(latticedefects, sclass, trials):
    def processor(ts):
        ts_new = ts.reduce_ends(checkld=latticedefects, _classes=tuple(sclass))
        return ts_new

    return processor


@alhambra.command(short_help="attempt to reduce the # of tiles")
@click.option("--latticedefects/--no-latticedefects", "-l/-L", default=True)
@click.option("--sclass", "-s", multiple=True, default=["2GO"])
@click.option("--rotate/--no-rotate", "-r/-R", default=True)
def reduce_tiles(latticedefects, sclass, rotate):
    def processor(ts):
        ts_new = ts.reduce_tiles(
            checkld=latticedefects, _classes=tuple(sclass), rotation=rotate
        )
        return ts_new

    return processor


@alhambra.command(short_help="optimize sticky end assignment")
def reorder_ends():
    pass


@alhambra.command(short_help="create strand sequences")
def create_strand_sequences():
    pass


@alhambra.command(short_help="run xgrow")
def run_xgrow():
    pass
