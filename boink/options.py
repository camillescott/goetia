import click
from boink import types
from boink import libboink

STORAGE_TYPES    = [t[0].__name__ for t in types.storage_types]
DEFAULT_ARGS = {t[0].__name__ : t[1] for t in types.storage_types}



@click.group(chain=True)
@click.option('--test', default='test')
def cli(test):
    click.echo(test)


def retrieve_storage(ctx, param, value):
    storage_t = getattr(libboink.storage, value)
    return storage_t


@cli.command()
@click.pass_context()
@cli.option('--store-type', type=click.Choice(STORAGE_TYPES), callback=retrieve_storage)
@cli.option('--store-args', nargs=-1)
def storage(ctx, store_type, store_args):
    print(store_type, store_args)


@cli.command()
@click.pass_context()
def stream(ctx):
    print('stream', ctx.store_type)
