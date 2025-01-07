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
    platforms=["Linux :: x86_64", "Windows :: x86_64", "Mac OS-X"],
    author="PenRed Contributors",
    author_email='vicent.gimenez.alventosa@gmail.com',
    url='https://github.com/penred/penred',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
    ],
    package_data={
        'pyPenred': ['*.so', '*.dll', '*.pyd', 'simulation*']
    },
    install_requires=[
        'pyyaml',
    ],
    python_requires='>=3.8',
    zip_safe=False,
)
