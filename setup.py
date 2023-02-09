import sys
from setuptools import setup, find_packages

# TODO: identify the purpose of this code block
# # update Git hash
# if len(sys.argv) > 1:
#     if sys.argv[1] == 'bdist_wheel':
#         try:
#             import git
#             repo = git.Repo(search_parent_directories=True)
#             sha = repo.head.object.hexsha
#             open('aeolis/GITVERSION', 'w').write(sha)
#         except ImportError:
#             print('*' * 70)
#             print('WARNING: Cannot update Git hash, because package "git" is not')
#             print('         installed. Continue packaging...')
#             print('*' * 70)
            

setup(
    name='AeoLiS',
    version=open('aeolis/VERSION').read().strip(),
    author='Sierd de Vries',
    author_email='sierd.devries@tudelft.nl',
    url='http://aeolis.readthedocs.io/',
    license='GNU GPLv3',
    description='A process-based model for simulating supply-limited aeolian sediment transport',
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords=['aeolian sediment transport coastal model deltares tudelft'],
    packages=find_packages(exclude=['docs', 'examples', 'tests', 'notebooks']),
    install_requires=[
        'docopt==0.6.1',
        'bmi-python',
        'netCDF4',
        'scipy',
        'numpy<1.24,>=1.18',
        'matplotlib',
        'numba'
    ],
    python_requires='>=3.7, <4',
    tests_require=[
        'pytest'
    ],
    test_suite='nose.collector',
    entry_points={'console_scripts': [
        'aeolis = aeolis.console:aeolis',
        'aeolis-wind = aeolis.console:wind',
    ]},
    include_package_data=True,
)
