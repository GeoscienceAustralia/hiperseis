import logging
import sys
import traceback
import warnings


def configure(verbosity):
    log = logging.getLogger("")
    log.setLevel(verbosity)
    ch = logging.StreamHandler()
    formatter = ElapsedFormatter()
    ch.setFormatter(formatter)
    log.addHandler(ch)


class ElapsedFormatter:

    def format(self, record):
        lvl = record.levelname
        name = record.name
        t = int(round(record.relativeCreated/1000.0))
        msg = record.getMessage()
        logstr = "+{}s {}:{} {}".format(t, name, lvl, msg)
        return logstr


def warn_with_traceback(message, category, filename, lineno, line=None):
    """
    copied from:
    http://stackoverflow.com/questions/22373927/get-traceback-of-warnings
    """
    traceback.print_stack()
    log = sys.stderr
    log.write(warnings.formatwarning(
        message, category, filename, lineno, line))
