from setuptools import setup, find_packages

setup(
    name="fastafunk",
    version="0.0.0",
    packages=find_packages(),
    url="https://github.com/cov-ert/fastafunk",
    license="MIT",
    entry_points={"console_scripts": ["fastafunk = fastafunk.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=[
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)