from setuptools import setup, find_packages

setup(
    name='sppran',
    version='0.0.1',
    description='A bioinformatic pipeline to analyze spatial proteomics samples using cell types defined by presence or abscent of protein markers',
    author='Sergio Zamora-Erazo',
    packages=find_packages(),
    install_requires=[
        'anndata>=0.10.8',
        'matplotlib>=3.9.1',
        'numpy>=2.0.0',
        'pandas>=2.2.2',
        'pyyaml>=6.0.1',
        'scanpy>=1.10.2',
        'scipy>=1.14.0',
        'seaborn>=0.13.2',
        'squidpy>=1.5.0',
        "statsmodels>=0.14.2"
    ],
    entry_points={
        'console_scripts': [
            'sppran=spatial_proteomics.core:main',  # optional CLI
        ],
    },
)