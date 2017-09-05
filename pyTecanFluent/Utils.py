from __future__ import print_function

# import
import os
import sys
import string
import logging
import itertools
import numpy as np

# functions
def file_written(file_name):
    """Status on writing file
    file_name: string
    """
    print('File written: {}'.format(file_name), file=sys.stderr)


def make_range(x, set_zero_index=False):
    """Making a range from a string of comma-delimited and hyphenated number strings.
    If x = 'all', None returned
    x : string
    set_zero_index : change from 1-index to 0-index
    Example: 1,2,3,4
    Example: 1,2,5-6
    """
    x = str(x)
    if x.lower() == 'all':
        return None
    y = str(x).split(',')

    z = []
    for i in y:
        j = [int(x) for x in i.split('-')]
        if len(j) > 1:
            z = z + list(range(j[0], j[1]+1))
        else:
            z = z + j

    if set_zero_index:
        z = [x-1 for x in z]
        if min(z) < 0:
            raise ValueError('Negative values from setting zero index')

    return z


def check_gwl(gwl_file):
    """Checking that gwl in correct format
    gwl_file : input file name or file handle
    """
    # file read
    try:
        gwl_file = open(gwl_file, 'r')
    except ValueError:
        pass
    # checks
    for line in gwl_file:
        line = line.rstrip()
        if line is '':
            msg = 'Empty lines not allowed in gwl files'
            raise ValueError(msg)
        if not line[0] in ('A', 'D', 'R', 'W', 'F', 'C', 'B'):
            msg = '"{}" not a valid command ID'
            raise ValueError(msg.format(line[0]))
    # file close
    gwl_file.close()


def to_win(file_name, suffix='_win'):
    """Create a copy of a file but with windows line breakds
    file_name : str, name of file
    suffix : added to file name of copy
    Returns : name of new file
    """
    x = os.path.splitext(file_name)
    out_file = x[0] + suffix + x[1]
    with open(file_name) as inFH, open(out_file, 'w') as outFH:
        for line in inFH:
            line = line.replace('\n', '\r\n')
            outFH.write(line)
    return out_file


def backup_file(f):
    """
    Back up a file, old_file will be renamed to #old_file.n#, where n is a
    number incremented each time a backup takes place
    """
    if os.path.exists(f):
        dirname = os.path.dirname(f)
        basename = os.path.basename(f)
        count = 1
        rn_to = os.path.join(
            dirname, '#' + basename + '.{0}#'.format(count))
        while os.path.exists(rn_to):
            count += 1
            rn_to = os.path.join(
                dirname, '#' + basename + '.{0}#'.format(count))
        logging.info("Backing up {0} to {1}".format(f, rn_to))
        os.rename(f, rn_to)
        return rn_to
    else:
        logging.warning('{0} doesn\'t exist'.format(f))


def position2well(position, wells=96, just_row=False, just_col=False):
    """Convert position to well
    Note: assuming column-wise ordering
    """
    # making plate index
    wells = int(wells)
    if wells == 96:
        nrows = 8
        ncols = 12
    elif wells == 384:
        nrows = 16
        ncols = 24
    else:
        raise ValueError('Number of wells ({}) not recognized'.format(wells))
    
    rows = list(string.ascii_uppercase[:nrows])
    cols = [x + 1 for x in range(ncols)]
    # position : [row, col]
    pos_idx = {i+1:x for i,x in enumerate(itertools.product(rows, cols))}
    # getting row-col
    try:
        row,col = pos_idx[position]
    except KeyError:
        msg = 'Cannot find well for position: {}'
        raise KeyError(msg.format(position))
    if just_row == True:
        return row
    elif just_col == True:
        return col
    else:
        return [row,col]
        
# main
if __name__ == '__main__':
    pass
