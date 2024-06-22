# setup.py
import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='SeaFreeze',
    version='1.0.0',
    author='Baptiste Journaux',
    author_email='bjournau@uw.edu',
    description='Thermodynamic properties of the phases of H2O and NaCl (aq)',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Bjournaux/SeaFreeze',
    download_url='https://pypi.org/project/SeaFreeze/',
    packages=['seafreeze', 'mlbspline', 'lbftd'],
    package_dir={'seafreeze': 'seafreeze', 'mlbspline': 'mlbspline', 'lbftd': 'lbftd'},
    classifiers=[
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent'
    ],
    include_package_data=True,
    install_requires=[
        'h5py >= 3.10',
        'hdf5storage >= 0.1.19',
        'numpy >= 1.26.3',
        'scipy >= 1.12',
        'psutil >= 5.9.8'
    ]
)
