import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="medipy",
    version="0.1.0",
    author="Nicolas Zabler",
    # author_email="neutral",
    description="Signal Analysis Package for Medical Data Analysis written in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zabler/medipy",
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn"
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
