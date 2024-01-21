from setuptools import setup, find_packages
import os
import shutil

def extract_Prior():
    if not os.path.exists('Prior'):
        os.makedirs('Prior')
        shutil.unpack_archive('Prior.tar.gz', '.', 'gztar')
        shutil.move('Prior', 'scPECA/Prior')

def custom_install():
    extract_Prior()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='scPECA',
    version='1.0',
    author='Jiahao Zhang',
    author_email='zhangjiahao@amss.ac.cn',
    description='PECA2 gene regulatory network construction for single-cell data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/zhangjiahao1234/scPECA',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # List of your package dependencies
    ],
    data_files=[('scPECA', ['Prior.tar.gz']), ('scPECA', ['Cones/*'])],
    cmdclass={'install': custom_install}
)