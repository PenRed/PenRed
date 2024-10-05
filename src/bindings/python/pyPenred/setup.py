from setuptools import setup, find_packages
import os

# Inclou el fitxer .so en el paquet
setup(
    name='pyPenred',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    description='Python interface for penRed framework',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    platforms=["linux-x86_64"],
    author='Vicent Gimenez Alventosa',
    author_email='vicent.gimenez.alventosa@gmail.com',
    url='https://github.com/penred/penred',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
    ],
    package_data={
        'pyPenred': ['*.so', '*.dll']
    },
    install_requires=[
        'pyyaml',
    ],
    zip_safe=False,
)
