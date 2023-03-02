import setuptools

setuptools.setup(
    name="medipy",
    version="0.0.6",
    author="Nicolas Zabler",
    author_email="neutralecho22@icloud.com",
    description="Code Collection for Medical Data Processing and Analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/zabler/medipy",
    packages=setuptools.find_packages('src', exclude=["*tests.*", "*tests", "*example.*"]),
    package_dir={'': 'src'},
    install_requires=[
        'numpy>=1.21.4',
        'scipy>=1.7.3',
        'matplotlib>=3.5.0',
    ],
    keywords=["medical", "signals", "analysis", "processing"],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
