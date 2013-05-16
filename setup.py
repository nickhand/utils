from setuptools import setup, find_packages
import os

setup(
    name='utils',
    version='1.1',
    author='Nick Hand',
    author_email='nicholas.adam.hand@gmail.com',
    packages=find_packages(),
    scripts=['bin/' + script for script in os.listdir('bin')],
    description='general utilities and gaussian processes regression module'
)