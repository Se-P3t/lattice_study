import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pylattice-Se_P3t",
    version="0.0.1",
    author="Se_P3t",
    author_email="sep3t.251124@gmail.com",
    description="lattice algorithms/attacks for CTF challenges",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Se-P3t/pylattice",
    packages=["pylattice", "pylattice.algorithms", "pylattice.attacks"],
    scripts=[],
    license='LICENSE',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.7',
    install_requires=["fpylll", "g6k"],
)