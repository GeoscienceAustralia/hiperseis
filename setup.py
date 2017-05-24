from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys


python_version = sys.version_info
__version__ = "0.0.1"

# numpy support for python3.3 not available for version > 1.10.1
if python_version.major == 3 and python_version.minor == 3:
    NUMPY_VERSION = 'numpy >= 1.9.2, <= 1.10.1'
else:
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

readme = open('README.md').read()

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
    packages=[''],
    # package_dir={'': ''},
    include_package_data=True,
    entry_points={
        'console_scripts': [
        ]
    },
    setup_requires=[NUMPY_VERSION],  # required due to netCDF4
    install_requires=[
        'Click >= 6.0',
        NUMPY_VERSION,
        'Cython >= 0.22.1',
        'mpi4py == 2.0.0',
        'scipy >= 0.15.1',
        'PyYAML >= 3.11',
        'matplotlib == 1.5.1',
        'pyproj >= 1.9.5',
        'Pillow >= 2.8.2',
        'joblib',
        'obspy >= 1.0.3',
        'h5py == 2.6.0',
        'pyasdf'
    ],
    extras_require={
        'dev': [
            'sphinx',
            'ghp-import',
            'sphinxcontrib-programoutput'
        ]
    },
    tests_require=[
        'pytest-cov',
        'coverage',
        'codecov',
        'tox',
        'pytest'  # pytest should be last
    ],
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
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
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
