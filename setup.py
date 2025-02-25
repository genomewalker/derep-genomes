from setuptools import setup, find_packages
import versioneer

requirements = [
    "pandas>=2.0.0",  # Updated for Python 3.11 compatibility
    "scipy>=1.10.0",  # Updated for Python 3.11 compatibility
    "networkx>=3.1",  # Updated for compatibility with newer Python versions
    "simple-slurm>=0.2.3",  # Maintained version
    "tqdm>=4.65.0",  # Latest stable version
    "PyYAML>=6.0",  # Latest stable version
    "leidenalg>=0.9.1",  # Updated for Python 3.11 compatibility
    "python-igraph>=0.10.4",  # Updated for Python 3.11 compatibility
    "biopython>=1.81",  # Latest stable version for Python 3.11
    "duckdb>=0.8.1",  # Added for the refactored database layer
    "numpy>=1.24.0",  # Added as it's a direct dependency
]

setup(
    name="derep-genomes",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A genome de-replication tool with fastANI, SKANI, Mash, and Dashing support",
    license="MIT",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/derep-genomes",
    packages=find_packages(),
    entry_points={"console_scripts": ["derepG=derep_genomes.__main__:main"]},
    install_requires=requirements,
    python_requires=">=3.8,<3.12",  # Specify Python version compatibility
    keywords="derep-genomes genome-dereplication fastani skani mash dashing",
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
    ],
)
