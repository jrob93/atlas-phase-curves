import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="atlasPhaseCurves",
    version="1.1.1",
    # author="James Robinson",
    # author_email="jrobinson72@qub.ac.uk",
    # description="A small example package",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    # project_urls={
    #     "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    # },
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
    # package_dir={"": "atlasPhaseCurves"},
    # packages=setuptools.find_packages(where="atlasPhaseCurves"),
    # package_dir={"": "atlas-phase-curves"},
    # packages=setuptools.find_packages("atlas-phase-curves"),
    # packages=setuptools.find_packages(),
    packages = ["atlas-phase-curves","atlas-phase-curves.calculate_phase","atlas-phase-curves.create_table"],
    # python_requires=">=3.9",
)
