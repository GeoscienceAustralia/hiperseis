import sys
from setuptools.command.test import test as TestCommand
from numpy.distutils.core import Extension, setup

python_version = sys.version_info
__version__ = "1.0.0"  # FZ-2020-07-02 tag this as version 1.0.0 after EFTF-1 completed and EFTF-X began.

# Read dependencies from requirements.in
with open('requirements.in', 'r') as f:
    install_reqs = f.readlines()
    # Remove comments and pip options (e.g. '--no-binary')
    install_reqs = [s.strip().split('#')[0].split('--')[0] for s in install_reqs]

print(install_reqs)

class PyTest(TestCommand, object):
    user_options = [('pytest-args=', 'a', "Arguments to pass to pytest")]

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

# The f2py fortran extension
# --------------------------
ellipcorr = Extension(
    name='ellipcorr',
    # add several modules files under the same extension
    sources=['ellip-corr/ellip/ellipcorr.f']
)

kennett_dist = Extension(
    name='kennett_dist',
    sources=['kennett-dist/kennett_dist.f']
)

setup(
    name='Passive-Seismic',
    version=__version__,
    description='Repository for development of software and '
                'metadata for passive seismic project',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Geoscience Australia Passive Seismic Team',
    author_email='',
    url='https://github.com/GeoscienceAustralia/Passive-Seismic',
    packages=['seismic',
              'seismic.traveltime',
              'seismic.receiver_fn',
              'seismic.gps_corrections',
              'seismic.amb_noise',
              'seismic.pick_harvester',
              'seismic.xcorqc',
              'seismic.ml_classifier',
              'seismic.ASDFdatabase',
              'seismic.inventory',
              'seismic.inversion',
              'seismic.synthetics',
              'PhasePApy'],
    package_dir={'passive-seismic': 'seismic'},
    include_package_data=True,
    dependency_links=[
    ],
    entry_points={
        'console_scripts': [
            'anulog = convert_logs.decode_datfile:anulog',
            'seismic = seismic.scripts.event:cli',
            'cluster = seismic.cluster.cluster:cli',
        ]
    },
    ext_modules=[ellipcorr, kennett_dist],
    # numpy preinstall required due to obspy
    # mpi4py  preinstall required due to h5py
    setup_requires=[
        'numpy',
        'mpi4py==3.0.0',
        'decorator>=4.1.0',
        'setuptools>=36.2.1'
    ],
    install_requires=install_reqs,
    extras_require={
        'dev': [
            'sphinx',
            'ghp-import',
            'sphinxcontrib-programoutput',
            'pytest-cov',
            'coverage==4.4.1',
            'codecov==2.0.9',
            'tox',
            'pytest>=3.1.0',
            'pytest-lazy-fixture>=0.4.0',
            'pytest-flake8>=0.8.1',
            'pytest-mock>=1.6.0',
            'pytest-cov==2.5.1',
            'pytest-regtest>=0.15.1',
            'flake8-docstrings>=1.1.0',
        ]
    },
    license="GNU GENERAL PUBLIC LICENSE v3",
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
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
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
