# medipy
Code Collection for Medical Data Processing and Analysis

[![Build Status](https://travis-ci.org/username/myproject.svg?branch=master)](https://travis-ci.org/username/myproject)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This Python package provides a collection of tools and utilities for processing and analyzing medical data. It includes modules to process data for different signal types like ECG, EEG, ACC, EMG, GYRO, etc. Every modul reflects an desribed model, algoeithm or concept published by a 

## Installation
You can install the package using pip:

```bash
pip install medipy
```

## Usage
To use the package, simply import the modules you need in your Python code: Here's an example of how to use a function by importing the complete package

```python
import medipy as mp
# Load ECG Data
r_peaks = mp.hamilton.function()
```

Here's another example of how to use a function by importing only the correspnding module

```python
from medipy import hamilton
# Load ECG Data
return = hamilton.function()
```

## License
This project is licensed under the MIT License - see the [License](LICENSE) file for details.



