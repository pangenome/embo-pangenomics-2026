# Practical Session

Guided exercises on genotyping from a pangenome using povu.


## Tooling

### Povu
1. Fetch the Code

Clone the repository and navigate into the project directory:

```
git clone https://github.com/urbanslug/povu.git
cd povu
```

2. Compile the project
Configures and compiles povu using cmake
```
cmake -H. -Bbuild && cmake --build build -- -j 3
```

3. Verify the Installation

Ensure that povu compiled successfully:
```
./bin/povu -h
```

### Bandage

- [Bandage](https://rrwick.github.io/Bandage/)

## Data

Fetch the data
```
git clone https://github.com/urbanslug/pang-datasets.git
```