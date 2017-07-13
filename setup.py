"""Credit: https://github.com/pypa/sampleproject"""
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.MD'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='optstoicpy',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.5.0',

    description='OptStoic python package',
    long_description=long_description,

    # The project's main homepage.
    url='http://www.maranasgroup.com/software.htm',

    # Author details
    author='Chiam Yu Ng',
    author_email='ngchiamyu@gmail.com',

    # Choose your license
    license='GNU GPLv3',

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
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',

    ],
    packages=find_packages(exclude= ['build',
                                     'data',
                                     'docs',
                                     'examples',
                                     'gams_mdf',
                                     'gams_optstoic',
                                     'mdf_ratio_script',
                                     'mdf_script',
                                     'obsolete',
                                     'optstoic_new_db',
                                     'proteincost',
                                     'res',
                                     'result',
                                     'result_v1',
                                     'tests',
                                     'util']),


    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[ 'pandas>=0.18.0',
                       #'scipy>=0.17.0',
                       #'numpy>=1.11.1',
                       'json>=2.0.9',
                       'graphviz>=0.4.8',
                       'PuLP>=1.6.1'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # extras_require={
    #     'dev': ['check-manifest'],
    #     'test': ['coverage'],
    # },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_dir={'optstoicpy': 'optstoicpy'},
    package_data={
        'optstoicpy': [ 'data/*.csv',
                        'data/*.json',
                        'data/optstoic_db_v3/*.txt',
                        'data/optstoic_db_v3/*.json',
                        'data/optstoic_db_v3/*.pkl'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('data', ['data/cofactors.csv',
    #                       'kegg_compound.json'])],

)