from setuptools import setup, find_packages

# Lê as dependências diretamente do requirements.txt
with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='atmos-bud',
    version='0.0.0',
    description='Program for analyzing the heat, vorticity and humidity budgets of limited regions on the atmosphere.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Danilo Couto de Souza',
    author_email='danilo.oceano@gmail.com',
    url='https://github.com/daniloceano/ATMOS-BUD',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=required,
    entry_points={
        'console_scripts': [
            'atmos-bud=atmos_bud:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    license='GPL-3.0',
    python_requires='>=3.10',
    license_files=('LICENSE',),
)
