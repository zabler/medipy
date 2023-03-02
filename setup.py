import setuptools

setuptools.setup(
    name="medipy",
    version="0.0.1",
    author="Nicolas Zabler",
    author_email="neutralecho22@icloud.com",
    description="Code Collection for Medical Data Processing and Analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/zabler/medipy",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy>=1.21.4',
        'scipy>=1.7.3',
    ],
    keywords=["medical", "signals", "algorithm"],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
