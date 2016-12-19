"""A bioinformatics pipeline for estimation of relative cell periods."""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import environ

environ['LC_ALL']='en_US.UTF-8'
environ['LANG']='en_US.UTF-8'

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md').read()

#here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
#with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

setup(
    name='menace',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.1.2',

    description='A metagenomics pipeline to estimate relative cell periods.',
    
    long_description=long_description,

    # The project's main homepage.
    url='https://www.github.com/SysBio/menace/',

    # Author details
    author='Daniel Hermansson',
    author_email='hedani@chalmers.se',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='systems biology research',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #py_modules=["ptr_pipeline.py"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    #setup_requires=[],
    install_requires=["numpy","scipy","pandas","biopython","matplotlib","xmltodict","configparser","lmfit","newick","Jinja2","docker"],
 #11 bx-python
 #-e git+https://github.com/PathoScope/PathoScope.git#egg=pathoscope],
    dependency_links = ['https://github.com/PathoScope/PathoScope/tarball'],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    #extras_require={
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #},

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'menace.test': ['comm00_1.fastq','comm00_2.fastq','searchStrings'],
        'menace.templates': ['jobscript','jobscript_local','project.conf'],
        'menace.bin': ['interp.pl','mainBuildbowtie2.sh','buildHelper.sh','mainBuildgem.sh','changeTID.sh'],
        'menace.extra': ['accLoc.csv','menace_run.ipynb'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'menace=menace.__main__:main',
            'fetch_seq=menace.bin.fetchSeq:main',
        ],
    }
)