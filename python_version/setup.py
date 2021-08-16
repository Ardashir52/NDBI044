from setuptools import setup

setup(
    name='bioinfobox',
    version=0.3,
    description='Bioinformatics toolbox',
    classifiers=[
        'Programming Language :: Python :: 3.8',
    ],
    url='https://github.com/Ardashir52/NDBI044/python_version',
    install_requires=[
        'biopython',
        'numpy'
    ],
    include_package_data=True,
    zip_safe=False,
    packages=[
        'bioinfobox',
    ],
    entry_points={
        'console_scripts': [
            'bioinfobox=bioinfobox.main:frontend'
        ]
    },
)
