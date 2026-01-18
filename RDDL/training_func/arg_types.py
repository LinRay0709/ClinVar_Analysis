import argparse
import os


def case_name(name):

    name = os.path.basename(os.path.abspath(name))

    if os.path.isdir(name):
        return name
    else:
        raise argparse.ArgumentTypeError('name of working directory %s does not exist.' % name)


def data_size(size):

    try:
        size = float(size)
    except ValueError as e:
        raise argparse.ArgumentTypeError(e)

    if 0 < size < 1:
        return size
    else:
        raise argparse.ArgumentTypeError('the percentage for data set splitting must be between 0.0 and 1.0.')
