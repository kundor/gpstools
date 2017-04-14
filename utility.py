"""Miscellaneous classes which are of use in gpsdata, particularly in rinex.py.

These are not very specific in usage, however, and could be useful anywhere.
"""
from contextlib import suppress, redirect_stdout, contextmanager
import subprocess
import os
import shutil
from textwrap import TextWrapper
from email.mime.text import MIMEText
from smtplib import SMTP
import traceback
import numpy as np
import config

def stacknames():
    istk = traceback.extract_stack()
    end = ''.join(f[2][0] for f in istk).rfind('<') + 1
    return ':'.join(f[2] for f in istk[end:-2])

def info(*args, logfile=None):
    if logfile is None:
        logfile = config.LOGFILE
    print(stacknames() + ':', *args, file=logfile, flush=True)

def debug(*args, **kwargs):
    if config.DEBUG:
        info(*args, **kwargs)

def vapr_email(subj, msg):
    """Email the string *msg* with subject *subj* to config.ERR_EMAIL from config.FROM_EMAIL."""
    msg = MIMEText(msg)
    msg['Subject'] = subj
    msg['From'] = config.FROM_EMAIL
    msg['To'] = config.ERR_EMAIL
    with SMTP(config.SMTP_HOST) as smtp:
        smtp.send_message(msg)

def email_exc(also_print=True):
    """Email a traceback of the exception currently being handled to config.ERR_EMAIL."""
    msg = traceback.format_exc()
    vapr_email('VAPR Error', msg)
    if also_print:
        info(msg)

def mode(arr):
    """Return the value occuring most often in numpy array arr.

    In case of a tie, returns the least value."""
    vals, ct = np.unique(arr, return_counts=True)
    return vals[np.argmax(ct)]

@contextmanager
def pushdir(ndir):
    """A context manager to change to the given directory, then change back when done."""
    old = os.getcwd()
    try:
        os.chdir(ndir)
    except FileNotFoundError:
        os.makedirs(ndir)
        os.chdir(ndir)
    try:
        yield
    finally:
        os.chdir(old)

@contextmanager
def stdouttofile(file):
    """A decorator to redirect stdout within a function to the given filename."""
    with open(file, 'a') as f, redirect_stdout(f):
        yield

@contextmanager
def fullprint(**kwargs):
    """A context manager to print full numpy arrays.

    with fullprint():
        print(a)
    with fullprint(linewidth=20, precision=3):
        print(a)
    """
    if 'threshold' not in kwargs:
        kwargs['threshold'] = np.inf
    oldopt = np.get_printoptions()
    np.set_printoptions(**kwargs)
    yield
    np.set_printoptions(**oldopt)

def static_vars(**kwargs):
    """A function decorator to set static variables as attributes on the function.

    Usage:
    @static_vars(counter=0)
    def foo():
        foo.counter += 1
        ...
    (User Claudiu, http://stackoverflow.com/a/279586/2132213)
    """
    def decorate(func):
        for k, val in kwargs.items():
            setattr(func, k, val)
        return func
    return decorate


def decompress(filename, move=False):
    """Decompress a (Lempel-Ziv) compress'd file.

    There seems to be no Python module to do this (the gzip module won't handle it,
    though the gzip program will), so we call an external process.
    These programs will only decompress if the filename ends with .Z;
    they remove the original file and output a file without the .Z.
    """
# However, given the -c flag, these programs will decompress to stdout, even if the filename
#  does not end with .Z.
    if not filename.endswith('.Z'):
        if move:
            defile = filename
            filename = filename + '.Z'
            os.rename(defile, filename)
        else:
            raise ValueError('Given filename ' + filename + ' does not end with .Z.')
    else:
        defile = filename[:-2]
    decompresscmds = [['compress', '-d', filename],
                      ['uncompress', filename],
                      ['gunzip', filename],
                      ['gzip', '-d', filename]]
    for cmd in decompresscmds:
        try:
            subprocess.run(cmd, check=True)
        except (OSError, subprocess.CalledProcessError):
            info("Command '", ' '.join(cmd), "' failed. Trying another...")
            continue
        if os.path.isfile(defile):
            return defile
        else:
            info("Command '", ' '.join(cmd), "' succeeded, but did not produce the output file?!")
    raise RuntimeError('Could not get an external program to decompress the file ' + filename)

def subdent(lead, obj, width):
    """Return lead + str(obj), with lines after the first indented to len(lead)."""
    if isinstance(obj, np.ndarray):
        obst = np.array2string(obj, max_line_width=width, prefix=lead)
    else:
        iin = ' '*len(lead)
        wrp = TextWrapper(width=width, expand_tabs=False, replace_whitespace=False,
                          initial_indent=iin, subsequent_indent=iin+'…')
        obst = '\n'.join(wrp.fill(txt) for txt in str(obj).splitlines()).lstrip()
    return lead + obst

class FieldObj:
    """Just an object you can set fields on."""
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        width = shutil.get_terminal_size().columns
        with fullprint(precision=3, threshold=20, infstr='∞'):
            return '\n'.join(subdent(k + ': ', v, width) for k, v in self.__dict__.items())

    def __repr__(self):
        pref = 'FieldObj('
        intr = ',\n    '
        return pref + intr.join('{}={!r}'.format(k, v) for k, v in self.__dict__.items()) + ')'

class ProfileThis:
    """A context manager to profile the contained code.

    Use:
    with ProfileThis() as pr:
        <code>
    p = pstats.Stats(pr)
    p.strip_dirs().sort_stats('cumulative').print_stats(20)
    """
    def __init__(self):
        import cProfile
        self.pr = cProfile.Profile()
    def __enter__(self):
        self.pr.enable()
        return self.pr
    def __exit__(self, *args):
        self.pr.disable()

class fileread(object):
    """A line-counting file wrapper.

    Wrap "sufficiently file-like objects" (ie those with readline())
    in an iterable which counts line numbers, strips newlines, and raises
    StopIteration at EOF.
    """
    def __new__(cls, file):
        """Create a fileread object.

        Input can be filename string, file descriptor number, or any object
        with `readline'.
        """
        if isinstance(file, fileread):
            file.reset()
            return file
        fr = object.__new__(cls)
        if isinstance(file, str):
            fr.fid = open(file)
            fr.name = file
        elif isinstance(file, int):
            fr.fid = os.fdopen(file)
            fr.name = "FD: " + str(file)
        elif hasattr(file, 'readline'):
            fr.fid = file
            if hasattr(file, 'name'):
                fr.name = file.name
            elif hasattr(file, 'url'):
                fr.name = file.url
        else:
            raise ValueError("Input of type " + str(type(file)) +
                             " is not supported.")
        fr.reset()
        return fr

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def next(self):
        """Return the next line, also incrementing `lineno'."""
        line = self.fid.readline()
        if not line:
            raise StopIteration()
        self.lineno += 1
        return line.rstrip('\r\n')

    __next__ = next

    def readline(self):
        """A synonym for next() which doesn't strip newlines or raise StopIteration."""
        line = self.fid.readline()
        if line:
            self.lineno += 1
        return line

    def __iter__(self):
        return self

    def reset(self):
        """Go back to the beginning if possible. Set lineno to 0 regardless."""
        if hasattr(self.fid, 'seek'):
            with suppress(OSError):
                self.fid.seek(0)
        self.lineno = 0

    def close(self):
        """Close the file.  A closed file cannot be used for further I/O."""
        if hasattr(self.fid, 'fileno') and self.fid.fileno() < 3:
            # Closing stdin, stdout, stderr can be bad
            return
        if hasattr(self.fid, 'close'):
            with suppress(OSError, EOFError):
                self.fid.close()
        elif hasattr(self.fid, 'quit'):
            with suppress(OSError, EOFError):
                self.fid.quit()


