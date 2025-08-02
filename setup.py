from setuptools import setup, find_packages

setup(
    name="iobrpy",
    version="0.1.2",
    packages=find_packages(),
    # packages=find_packages(),  # 自动发现并包含所有包
    include_package_data=True,
    package_data={
        'iobrpy': ['resources/*'],  # 这样包含resources目录下的文件
    },
    entry_points={
        'console_scripts': [
            'iobrpy=iobrpy.main:main',   # 创建一个命令行工具链接到 main.py 中的 main 函数
        ],
    },
    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'matplotlib',
        'scikit-learn',
    ],
    author="Dongqiang Zeng",
    author_email='interlaken@smu.edu.cn',
    license='MIT License',
    description="Intratumoral Microbiome Finder Tool",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/IOBR/IOBRpy",
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
