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
        'bmi',
        'scipy',
        'numpy',
        'docopt',
    ],
    dependency_links=[
        'git+https://github.com/openearth/bmi-python.git#egg=bmi',
    ],
    tests_require=[
        'nose'
    ],
    test_suite='nose.collector',
    entry_points={'console_scripts': [
        'aeolis = aeolis.console:aeolis',
        'aeolis-wind = aeolis.console:wind',
    ]},
)
