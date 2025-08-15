
Getting Started
===============

This section provides a step-by-step guide on how to set up a development environment and get ATMOS-BUD up and running on your system.

Cloning the Repository
----------------------

Before setting up the environment, the first step is to clone the ATMOS-BUD repository from GitHub. Use the following `git` command to clone the repository:

.. code-block:: bash

    # Clone the ATMOS-BUD repository
    git clone https://github.com/daniloceano/ATMOS-BUD.git

    # Navigate to the repository directory
    cd ATMOS-BUD

Creating a Conda Environment
----------------------------

ATMOS-BUD requires a dedicated Python environment to manage dependencies. We recommend using Conda, a popular package, dependency, and environment management tool. To create a new Conda environment and activate it, follow these steps after you've cloned the repository and navigated to the directory:

.. code-block:: bash

    # Create a new Conda environment named 'atmosbud'
    conda create --name atmosbud python=3.10

    # Activate the environment
    conda activate atmosbud

Installing Dependencies
-----------------------

Once the Conda environment is set up and activated, install the required dependencies using the `requirements.txt` file provided in the repository. Run the following command to install all dependencies via `pip`:

.. code-block:: bash

    # Ensure pip, wheel, and setuptools are updated
    pip install --upgrade pip wheel setuptools

    # Install the required packages
    pip install -r requirements.txt

You are now ready to begin using ATMOS-BUD for your atmospheric science calculations and analysis.