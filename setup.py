from setuptools import setup, find_packages

setup(
    name='sppran',
    version='0.2.5',
    description='A bioinformatic pipeline to analyze spatial proteomics samples using cell types defined by presence or abscent of protein markers',
    author='Sergio Zamora-Erazo',
    packages=find_packages(),
    install_requires=[
        'anndata>=0.12.10',
        'matplotlib>=3.10.8',
        'numpy>=2.4.2',
        'pandas>=2.3.3',
        'pyarrow>=24.0.0',
        'pyyaml>=6.0.3',
        'scanpy>=1.12',
        'scipy>=1.17.1',
        'seaborn>=0.13.2',
        'squidpy>=1.8.1',
        "statsmodels>=0.14.6"
    ],
    entry_points={
        'console_scripts': [
            'sppran=spatial_proteomics.core:main',  # optional CLI
        ],
    },
)