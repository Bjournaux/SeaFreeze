import setuptools

setuptools.setup(
    name='SeaFreeze',
    version='0.1a1',
    author='pennythewho',
    author_email='who@pennythewho.com',
    description='thermodynamic properties of the phases of H₂O',
    # currently "hidden" in the seafreeze branch
    url='https://github.com/jmichaelb/PublicationArchive/tree/seafreeze/SeaFreeze',
    packages=['seafreeze'],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ]
)