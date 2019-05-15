#!/usr/bin/env python
"""
Helper functions to convert FDSN station XML files to Seiscomp3 SC3ML format.

Can be used as a standalone tool as well:
    `fdsnxml_convert.py src_path dst_path`
"""

import os
import sys
import subprocess
# import tempfile

import click
# from seismic.inventory.response import ResponseFactory
# from obspy import read_inventory

sc3_converter_app = "fdsnxml2inv"
sc3_converter_options = ("--quiet", "--formatted")


def toSc3ml(src_path, dst_path):
    """
    Convert file(s) in src_path from FDSN station XML to SC3ML and emit result(s) to dst_path.

    If src_path is a file, dst_path will be treated as a file. If dst_path already exists as a
    folder, an exception is raised.

    If src_path is a folder, dst_path will be treated as a folder. If dst_path already exists as
    a file, an exception is raised. The src_path directory hierarchy will be walked to find all
    .xml files, each of which will be converted to a mirrored relative path under dst_path.

    :param src_path: The source path from which to input XML file(s).
    :type src_path: str or pathlib.Path
    :param dst_path: The destination path into which converted sc3ml file(s) are output.
    :type dst_path: str or pathlib.Path
    :raises: OSError, FileNotFoundError, RuntimeError, FileExistsError
    """
    if not sc3_conversion_available():
        raise OSError(2, "No such application found on path", sc3_converter_app)

    if not os.path.exists(src_path):
        raise FileNotFoundError(src_path)

    # response = None
    # if response_fdsnxml is not None:
    #     rf = ResponseFactory()
    #     rf.CreateFromStationXML('resp', response_fdsnxml)
    #     response = rf.getResponse('resp')

    if os.path.isfile(src_path):
        # _fileToSc3ml(src_path, dst_path, response)
        _fileToSc3ml(src_path, dst_path)
        _reportConversion([src_path], [])
    elif os.path.isdir(src_path):
        # _reportConversion(*_folderToSc3ml(src_path, dst_path, response))
        _reportConversion(*_folderToSc3ml(src_path, dst_path))
    else:
        raise RuntimeError("Unsupported file type for {}".format(src_path))


def _fileToSc3ml(src_file, dst_file):
    assert os.path.isfile(src_file)
    if os.path.exists(dst_file) and os.path.isdir(dst_file):
        raise FileExistsError("{} already exists as a folder".format(dst_file))

    # # Repair dates and inject responses if needed.
    # inv = read_inventory(src_file)
    # for n in inv.networks:
    #     for s in n.stations:
    #         if(not s.start_date): s.start_date = n.start_date
    #         if(not s.end_date): s.end_date = n.end_date
    #         for c in s.channels:
    #             if(not c.start_date): c.start_date = s.start_date
    #             if(not c.end_date): c.end_date = s.end_date
    #             if(response and not c.response): c.response = response

    #     fn = os.path.join(tempfile.gettempdir(), os.path.basename(src_file))
    #     inv.write(fn, format='STATIONXML')
    #     src_file = fn

    cmd = [sc3_converter_app] + list(sc3_converter_options) + [src_file, dst_file]
    # Convert using system call
    print(" ".join(cmd))
    if sys.version_info[0] < 3:
        subprocess.check_call(cmd)
    else:
        subprocess.check_call(cmd, timeout=3600)


def _makedirs(path, exist_ok=True):
    if (sys.version_info >= (3, 0)):
        os.makedirs(path, exist_ok=exist_ok)
    else:
        if not os.path.exists(path): os.makedirs(path)


def _folderToSc3ml(src_folder, dst_folder):
    assert os.path.isdir(src_folder)
    if os.path.exists(dst_folder) and os.path.isfile(dst_folder):
        raise FileExistsError("{} already exists as a file".format(dst_folder))

    _makedirs(dst_folder, exist_ok=True)

    success_files = []
    failed_files = []
    for root, _, files in os.walk(src_folder):
        for f in files:
            relpath = os.path.relpath(root, src_folder)
            src_file = os.path.join(root, f)
            if relpath != '.':
                dst_tree = os.path.join(dst_folder, relpath)
            else:
                dst_tree = dst_folder
            _makedirs(dst_tree, exist_ok=True)
            dst_file = os.path.join(dst_tree, f)
            try:
                # _fileToSc3ml(src_file, dst_file, response)
                _fileToSc3ml(src_file, dst_file)
                success_files.append(src_file)
            except (subprocess.CalledProcessError, OSError) as e:
                failed_files.append((src_file, str(e)))

    return success_files, failed_files


def sc3_conversion_available():
    if (sys.version_info >= (3, 0)):
        from shutil import which
        return which(sc3_converter_app) is not None
    else:
        return True # assume it exists


def _reportConversion(success_files_list, failed_files_list):
    num_success = len(success_files_list)
    num_failed = len(failed_files_list)
    if num_success > 0:
        print("Successfully converted {} files to sc3ml".format(num_success))
    if num_failed > 0:
        print("Following {} files failed conversion:".format(num_failed))
        for f, e in failed_files_list:
            print("  {} ({})".format(f, e))


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('src_path', type=click.Path('r'))
@click.argument('dst_path', type=str)
def main(src_path, dst_path):
    toSc3ml(src_path, dst_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
