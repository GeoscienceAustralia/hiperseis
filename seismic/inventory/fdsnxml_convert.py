#!/usr/bin/env python
"""
Helper functions to convert FDSN station XML files to Seiscomp3 SC3ML format.

Can be used as a standalone tool as well:
    `fdsnxml_convert.py src_path dst_path`
"""

import os
import sys
import subprocess

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
    if not _checkConverterAvailability():
        raise OSError(2, "No such application found on path", sc3_converter_app)

    if not os.path.exists(src_path):
        raise FileNotFoundError(src_path)

    if os.path.isfile(src_path):
        _fileToSc3ml(src_path, dst_path)
        _reportConversion([src_path], [])
    elif os.path.isdir(src_path):
        _reportConversion(*_folderToSc3ml(src_path, dst_path))
    else:
        raise RuntimeError("Unsupported file type for {}".format(src_path))


def _fileToSc3ml(src_file, dst_file):
    assert os.path.isfile(src_file)
    if os.path.exists(dst_file) and os.path.isdir(dst_file):
        raise FileExistsError("{} already exists as a folder".format(dst_file))

    cmd = [sc3_converter_app] + list(sc3_converter_options) + [src_file, dst_file]
    # Convert using system call
    print(" ".join(cmd))
    if sys.version_info[0] < 3:
        subprocess.check_call(cmd)
    else:
        subprocess.check_call(cmd, timeout=3600)


def _folderToSc3ml(src_folder, dst_folder):
    assert os.path.isdir(src_folder)
    if os.path.exists(dst_folder) and os.path.isfile(dst_folder):
        raise FileExistsError("{} already exists as a file".format(dst_folder))

    os.makedirs(dst_folder, exist_ok=True)

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
            os.makedirs(dst_tree, exist_ok=True)
            dst_file = os.path.join(dst_tree, f + ".sc3ml")
            try:
                _fileToSc3ml(src_file, dst_file)
                success_files.append(src_file)
            except subprocess.CalledProcessError:
                failed_files.append(src_file)

    return success_files, failed_files


def _checkConverterAvailability():
    from shutil import which
    return which(sc3_converter_app) is not None


def _reportConversion(success_files_list, failed_files_list):
    num_success = len(success_files_list)
    num_failed = len(failed_files_list)
    if num_success > 0:
        print("Successfully converted {} files to sc3ml".format(num_success))
    if num_failed > 0:
        print("Following {} files failed conversion:".format(num_failed))
        for f in failed_files_list:
            print("  " + f)


if __name__ == "__main__":
    src = sys.argv[1]
    dst = sys.argv[2]
    print("Command: convert {} to {}".format(src, dst))
    toSc3ml(src, dst)
