from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name="derep-genomes",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A simple genome de-replication tool with fastANI",
    license="GNUv3",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/derep-genomes",
    packages=["derep_genomes"],
    entry_points={"console_scripts": ["derep-genomes=derep_genomes.__main__:main"]},
    install_requires=requirements,
    keywords="derep-genomes",
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
)
