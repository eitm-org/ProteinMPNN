#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pip

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

## workaround derived from: https://github.com/pypa/pip/issues/7645#issuecomment-578210649
parsed_requirements = parse_requirements(
    'requirements/prod.txt',
    session='workaround'
)

parsed_test_requirements = parse_requirements(
    'requirements/test.txt',
    session='workaround'
)


requirements = [str(ir.requirement) for ir in parsed_requirements]
test_requirements = [str(tr.requirement) for tr in parsed_test_requirements]


setup(
    name='proteinmpnn',
    version='0.1.0',
    description="Python Boilerplate contains all the boilerplate you need to create a Python package.",
    long_description=readme + '\n\n' + history,
    author="Xingyao Chen",
    author_email='xchen@eitm.org',
    url='https://github.com/xingyaoc/proteinmpnn',
    packages=[
        'proteinmpnn',
        'proteinmpnn.helper_scripts',
    ],
    package_dir={'proteinmpnn':
                 'proteinmpnn',
                 'proteinmpnn.helper_scripts': 'proteinmpnn/helper_scripts',
                 },
    include_package_data=True,
    install_requires=requirements,
    license="ISCL",
    zip_safe=False,
    keywords='proteinmpnn',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
