import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="theovib",
    version="0.0.1",
    author="L. J. Duarte",
    author_email="leo.j.duarte49@gmail.com",
    description="Tools for molecular vibrations analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ljduarte/theovib",
    project_urls={
        "Bug Tracker": "https://github.com/ljduarte/theovib/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)