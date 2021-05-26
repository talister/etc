from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["astropy>=4.0", "toml", "synphot>=0.3"]

setup(
    name="etc",
    version="0.0.3",
    author="Tim Lister",
    author_email="tlister@lco.global",
    description="A package to model astronomical sites, telescopes and instrument and predict exposure times",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/talister/etc/",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
