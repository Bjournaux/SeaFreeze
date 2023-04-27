import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='SeaFreeze',
    version='0.9.3',
    author='Marshall J. Styczinski',
    author_email='marshall.j.styczinski@jpl.nasa.gov',
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
    package_data={'seafreeze': ['SeaFreeze_Gibbs.mat', 'LICENSE.txt']},
    include_package_data=True,
    install_requires=['numpy', 'scipy', 'psutil']
)
