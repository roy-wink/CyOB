from setuptools import setup
# from setuptools import find_packages

setup(
    name='cyob',
    version='0.1.2',
    description='Python package for easily calculating and visualizing molecular orbitals.',
    url='https://github.com/roy-wink/CyOB',
    author='Roy Wink',
    author_email='r.wink@student.tue.nl',
    license='GPL3',
    packages=['cyob'],
    python_requires='>=3.5',
    install_requires=['numpy', 'matplotlib', 'setuptools']
)