# -*- coding: utf-8 -*-

import os
import os.path
from setuptools import find_packages
from setuptools import setup

name = 'cardio-clean'
version = '0.0.1'


def find_requires():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open('{0}/requirements.txt'.format(dir_path), 'r') as reqs:
        requirements = reqs.readlines()
    return requirements


if __name__ == "__main__":
    setup(
        name=name,
        version=version,
        description='Cardio signal enhancement',
        long_description="Scripts for ECG enhancement",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Programming Language :: Python'
        ],
        packages=find_packages(),
        install_requires=find_requires(),
        data_files=[],
        include_package_data=True,
        entry_points={
            'console_scripts': [
                'cardioclean = cardio_clean.cardioclean:main'
            ],
        },
    )
