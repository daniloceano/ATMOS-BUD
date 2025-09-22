from setuptools import setup, find_packages
import os

# Lê as dependências diretamente do requirements.txt
def read_requirements():
    with open('requirements.txt') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Lê o README
def read_readme():
    if os.path.exists('README.md'):
        with open('README.md', encoding='utf-8') as f:
            return f.read()
    return ''

setup(
    name='atmos-bud',
    version='0.1.6',
    description='Program for analyzing the heat, vorticity and humidity budgets of limited regions on the atmosphere.',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    author='Danilo Couto de Souza',
    author_email='danilo.oceano@gmail.com',
    url='https://github.com/daniloceano/ATMOS-BUD',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=read_requirements(),
    entry_points={
        'console_scripts': [
            'atmos-bud=atmos_bud.main:main',  # Verificar se este path está correto
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Operating System :: OS Independent',
    ],
    license='GPL-3.0',
    python_requires='>=3.10',
    license_files=('LICENSE',),
    keywords='atmospheric science, meteorology, budget analysis, ERA5',
    project_urls={
        'Documentation': 'https://github.com/daniloceano/ATMOS-BUD/docs',
        'Source': 'https://github.com/daniloceano/ATMOS-BUD',
        'Tracker': 'https://github.com/daniloceano/ATMOS-BUD/issues',
    },
)