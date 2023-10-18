import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='SeaFreeze',
    version='0.9.6',
    author='Baptiste Journaux',
    author_email='bjournau@uw.edu',
    description='Thermodynamic properties of the phases of H2O',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Bjournaux/SeaFreeze',
    packages=['seafreeze', 'mlbspline', 'lbftd'],
    package_dir={'seafreeze': 'seafreeze', 'mlbspline': 'mlbspline', 'lbftd': 'lbftd'},
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent'
    ],
    include_package_data=True,
    install_requires=[
        'numpy >= 1.24.2',
        'scipy >= 1.10.1',
        'psutil >= 5.9.4'
    ]
)
