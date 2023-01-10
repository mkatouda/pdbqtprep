from setuptools import setup

setup(
    name="pdbqtprep",
    version="0.0.1",
    install_requires=[
        'pyyaml', 'numpy', 'pdbfixer', 'openmm >= 7.1'
    ],
    entry_points={
        'console_scripts': [
            'pdbqtprep=pdbqtprep.pdbqtprep:main',
        ],
    },
    author="Michio Katouda",
    author_email="katouda@rist.or.jp",
    description="Autodock Vina protein PDBQT file prepare tool",
    url="https://github.com/mkatouda/pdbqtprep",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.7',
)
