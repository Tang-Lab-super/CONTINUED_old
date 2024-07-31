from setuptools import Command, find_packages, setup

__lib_name__ = "continued"
__lib_version__ = "1.0.0"
__description__ = "Cluster, integration and annotation of In situ metabonomics data"
__url__ = "https://github.com/Tang-Lab-super/CONTINUED"
__author__ = "Yuchen Yuan"
__author_email__ = "yuanych9@mail2.sysu.edu.cn"
__license__ = "MIT"
__keywords__ = ["In situ metabonomics"]
__requires__ = ["requests",]

with open("README.rst", "r", encoding="utf-8") as f:
    __long_description__ = f.read()

setup(
    name = __lib_name__,
    version = __lib_version__,
    description = __description__,
    url = __url__,
    author = __author__,
    author_email = __author_email__,
    license = __license__,
    packages = ['continued'],
    install_requires = __requires__,
    zip_safe = False,
    include_package_data = True,
    long_description = __long_description__
)