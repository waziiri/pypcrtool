from setuptools import setup, find_packages

setup(
    name="pypcrtool",
    version="1.1.1",
    description="A Python package for In silico PCR and primer verification.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Ibrahim Zubairu Waziri, Mustapha Ibrahim Usman, Zainab Ali Dandalma",
    author_email="biotechizwaziri@gmail.com, musteengumel@polac.edu.ng, zainabalidandalma92@gmail.com",
    url="https://github.com/waziiri/pypcrtool",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        "in silico PCR", 
        "primer specificity", 
        "gel electrophoresis", 
        "DNA amplification", 
        "bioinformatics",
        "computational biology"
    ],
    install_requires=[
        "numpy",
        "matplotlib"
    ],
    license="MIT",
    project_urls={
        "Source": "https://github.com/waziiri/pypcrtool",
    },
    python_requires='>=3.6',
)

