from setuptools import find_packages, setup

setup(
    name="pysurfacefit",
    version="0.1",
    author="Roman Kochanov",
    author_email="",
    description="Multidimensional data fitting tool in Python",
    #url="",
    python_requires=">=3.7",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "numba",
        "tabulate",
        "scipy",
        "matplotlib",
        "jeanny3",
    ],
    entry_points = {
        'console_scripts': ['pysurfacefit=pysurfacefit.command_line:main']
    }
)
