# Getting Started

## Prerequisites
- Python 3.10+
- CMake toolchain for C++ builds
- Optional: virtual environment for Python tooling

## Clone and build
```bash
git clone https://github.com/pnnl/qasmtrans.git
cd qasmtrans

python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

mkdir -p build
cmake -S . -B build
cmake --build build
```

## Run the transpiler
```bash
./build/qasmtrans -i data/test_benchmark/bv10.qasm -m ibmq -c data/devices/ibmq_toronto.json -v 1
```

## Python bindings
```bash
source venv/bin/activate
pip install -r requirements.txt
cmake -S . -B build
cmake --build build --target qasmtrans_core
```
Or install the Python package into the environment:
```bash
pip install .
```
