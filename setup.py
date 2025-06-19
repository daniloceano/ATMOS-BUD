from setuptools import setup, find_packages

setup(
    name='atmos-bud',  
    version='0.0.0', 
    description='Program for analyzing the heat, vorticity and humidity budgets of limited regions on the atmosphere.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Seu Nome',
    author_email='seuemail@exemplo.com',
    url='https://github.com/daniloceano/ATMOS-BUD',
    packages=find_packages(where='src'),
    package_dir={'': 'src'}, 
    install_requires=[
        'pandas==2.2.3',
        'numpy==1.26.4',
        'xarray==2025.4.0',
        'MetPy==1.6.1',
        'cdsapi>=0.6',
        'Cartopy==0.22.0',
        'dask==2024.2.0',
        'scipy==1.12.0'

    ],
    entry_points={
        'console_scripts': [
            'atmos-bud=atmos_bud:main',  # Ajuste de acordo com seu ponto de entrada
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
