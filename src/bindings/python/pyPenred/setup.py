from setuptools import setup, find_packages, Extension
import os

# Read the version from an environment variable or use the default
VERSION = os.getenv("PENRED_VERSION", "1.13.0")

ext_modules = [
    Extension(
        name="pyPenred.simulation",
        sources=[],  # Empty because the .so file is pre-built
    ),
]

setup(
    name='pyPenred',
    version=VERSION,
    packages=find_packages(),
    ext_modules=ext_modules,
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
