# setup.py
from setuptools import setup, find_packages

setup(
    name="dftoy",
    version="1.1.3",
    author="Julian Cogua, Juan Esteban Neira",  # varios autores separados por coma
    author_email="ojcoguaa@udistrital.edu.co, jueneirad@udistrital.edu.co",  # separados por coma
    description="Simulación DFT 1D con GUI",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/XxyamzxX/DFTOY",  # Cambia a tu repo si lo tienes
    packages=find_packages(),
    include_package_data=True,  # necesario para package_data
    package_data={
        "dftoy": ["Manual_Usuario/*.pdf"],  # incluye todos los PDFs de la carpeta Manual_Usuario
    },
    install_requires=[
        "numpy",
        "matplotlib",
    ],
    entry_points={
        "console_scripts": [
            "dftoy-gui=dftoy.gui:main",  # Ejecuta la GUI
            "dftoy=dftoy.main:main",     # Ejecuta el main.py
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.11",
)

