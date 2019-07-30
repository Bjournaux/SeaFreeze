import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='SeaFreeze',
    version='0.8.1a',
    author='pennythewho',
    author_email='who@pennythewho.com',
    description='thermodynamic properties of the phases of H₂O',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/Bjournaux/SeaFreeze',
    install_requires=['uw-highP-geophysics-tools>=0.8'],
    packages=['seafreeze'],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    package_data={'seafreeze': ['SeaFreeze_Gibbs.mat'], '': ['LICENSE.txt']},
    include_package_data=True
)