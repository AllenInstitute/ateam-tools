from setuptools import setup, find_packages
setup(
    name='ateam',
    version='0.0.1',
    author='Tom Chartrand',
    author_email='tom dot chartrand at alleninstitute dot org',
    description='A-team python tools',
    install_requires=[
        'pandas',
        'numpy',
        'bmtk[bionet]'
    ],
    packages=find_packages(),
    scripts=['bin/run_hpc'],
    platforms='any'
)
