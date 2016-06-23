from setuptools import setup, find_packages

setup(
    name='AeoLiS',
    version='0.0',
    author='Bas Hoonhout',
    author_email='b.m.hoonhout@tudelft.nl',
    packages=find_packages(),
    description='A process-based model for simulating supply-limited aeolian sediment transport',
    long_description=open('README.txt').read(),
    install_requires=[
        'numpy',
        'scipy',
        'docopt',
        'bmi',
        'netCDF4',
    ],
    setup_requires=[
        'sphinx',
        'sphinx_rtd_theme'
    ],
    tests_require=[
        'nose'
    ],
    test_suite='nose.collector',
    entry_points={'console_scripts': [
        'aeolis = aeolis.cmd:aeolis',
        'aeolis-wind = aeolis.cmd:wind',
    ]},
    data_files=[
        ('example', [
            "examples/",
        ])
    ],
)
