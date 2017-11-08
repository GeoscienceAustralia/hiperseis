from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys


python_version = sys.version_info
__version__ = "0.0.1"

NUMPY_VERSION = 'numpy >= 1.9.2'


class PyTest(TestCommand, object):

    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        super(PyTest, self).initialize_options()
        self.pytest_args = []

    def finalize_options(self):
        super(PyTest, self).finalize_options()
        self.test_suite = True
        self.test_args = []

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        exit(pytest.main(self.pytest_args))

readme = open('README.rst').read()

doclink = """
Documentation
-------------

The full documentation is at http://geoscienceaustralia.github.io/passive
-seismic
/."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='Passive-Seismic',
    version=__version__,
    description='Repository for development of software and '
                'metadata for passive seismic project',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Geoscience Australia Passive Seismic Team',
    author_email='',
    url='https://github.com/GeoscienceAustralia/Passive-Seismic',
    packages=['seismic', 'convert_logs', 'seismic.scripts', 'seismic.cluster'],
    package_dir={'passive-seismic': 'seismic'},
    include_package_data=True,
    dependency_links=[
        'https://github.com/matplotlib/basemap/archive/'
        'v1.1.0.zip#egg=basemap-1.1.0',
        # 'https://github.com/GeoscienceAustralia/PhasePApy/'
        # 'archive/master.zip#egg=phasepapy-1.1.1',
        # git+https://github.com/GeoscienceAustralia/PhasePApy.git@v1.1.1
        ],
    entry_points={
        'console_scripts': [
            'anulog = convert_logs.decode_datfile:anulog',
            'seismic = seismic.scripts.event:cli',
            'cluster = seismic.cluster.cluster:cli',
        ]
    },

    # numpy preinstall required due to obspy
    # mpi4py  preinstall required due to h5py
    setup_requires=[
        NUMPY_VERSION,
        'mpi4py==2.0.0',
        'decorator>=4.1.0',
        'setuptools>=36.2.1'
    ],
    install_requires=[
        'Click >= 6.0',
        NUMPY_VERSION,
        'Cython >= 0.22.1',
        'mpi4py == 2.0.0',
        'scipy >= 0.15.1',
        'PyYAML >= 3.11',
        'matplotlib >= 1.4.3',  # nci version with python=3.4
        'joblib',
        'obspy >= 1.0.3',  # 1.0.3 does not have sc3ml read functionality
        'h5py >= 2.6.0',
        'pyasdf',
        'pyqtgraph',
        # 'phasepapy == 1.1.1',
        'basemap == 1.1.0',
        'netCDF4 >= 1.3.0',
    ],
    extras_require={
        'dev': [
            'sphinx',
            'ghp-import',
            'sphinxcontrib-programoutput',
            'pytest-cov',
            'coverage == 4.4.1',
            'codecov == 2.0.9',
            'tox',
            'pytest >= 3.1.0',
            'pytest-flake8 >= 0.8.1',
            'pytest-mock >= 1.6.0',
            'pytest-cov == 2.5.1',
            'pytest-regtest >= 0.15.1',
            'flake8-docstrings >= 1.1.0',
        ]
    },
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='Passive Seismic',
    classifiers=[
        'Development Status :: 4 - Beta',
        "Operating System :: POSIX",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        # "Programming Language :: Python :: 3",
        # "Programming Language :: Python :: 3.3",
        # "Programming Language :: Python :: 3.4",
        # "Programming Language :: Python :: 3.5",
        # "Programming Language :: Python :: 3.6",
        # "Programming Language :: Python :: 3.7",
        # add additional supported python versions
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Information Analysis"
        # add more topics
    ],
    cmdclass={
        'test': PyTest,
    }
)
