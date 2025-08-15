Output Management Module
========================

The output management module handles saving and exporting analysis results in various formats.

Main Functions
--------------

.. automodule:: src.output_management
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Key Features
------------

* CSV output generation
* NetCDF file creation
* Result formatting
* Metadata preservation
* File organization
* Data compression

Output Formats
--------------

CSV Files
~~~~~~~~~

Budget results exported as comma-separated values:

* Time series data
* Statistical summaries
* Parameter tables
* Processing logs

NetCDF Files
~~~~~~~~~~~~

Structured output preserving:

* Original coordinate systems
* Variable attributes
* Processing history
* Dimensional relationships

Directory Structure
~~~~~~~~~~~~~~~~~~~

Organized output directory:

* Results by domain type
* Timestamped analysis runs
* Parameter-specific subdirectories
* Summary files

Usage Examples
--------------

.. code-block:: python

   from src.output_management import save_results, export_csv
   
   # Save complete results
   save_results(budget_results, output_dir='Results/')
   
   # Export specific data to CSV
   export_csv(heat_budget_data, 'heat_budget_results.csv')
