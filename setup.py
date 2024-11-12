from setuptools import setup, find_packages

setup(
    name="menthu",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "biopython==1.79",
        "pandas==1.3.5",
        "numpy==1.21.6",
        "scikit-learn==0.20.0",
        "scipy==1.7.3",
        "joblib>=0.13.2,<0.14.0",  # Python 3.7対応バージョン
        "packaging==21.3",
    ],
    entry_points={
        'console_scripts': [
            'menthu=menthu.menthu:main',
        ],
    },
    author="Ryo Niwa",
    description="MENTHU (Microhomology-mediated End joining Target HSting Utility)",
    python_requires=">=3.7, <3.8",
)