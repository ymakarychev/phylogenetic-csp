# Arity-4 computations for satisfiable phylogenetic CSPs

This repository contains the code and data for the arity-4 results in Yury Makarychev's *Approximation algorithms for satisfiable phylogenetic CSPs*: the approximation guarantees in Table 1 and Appendix C, and the predicate catalogue in Appendix D.

See [results.md](results.md) for the guarantees, satisfaction probabilities, and catalogue. The predicate IDs and pattern names agree with those in the paper.

Run the commands below from the repository directory.

## Regenerate the predicate classification

Compile the classifier with a C++20 compiler, such as `clang++` or `g++`:

```sh
clang++ -std=c++20 -O2 -o csp_classification csp_classification.cpp
```

Run it and save the output:

```sh
./csp_classification > arity4/classification.txt
```

The classifier enumerates all 32,768 predicates and groups them into 1,576 equivalence classes. Of these, 12 are tractable (including the constant predicates), and 40 are NP-hard with a canonical relaxation that is not always true. The saved output lists these 40 predicates and their relaxations.

## Verify the approximation guarantees

With Python 3.10 or later, create and activate a virtual environment and install the dependencies:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the verification script:

```sh
python verify_arity4.py \
    --classification arity4/classification.txt \
    --schemes arity4/rounding_schemes.py
```

The verifier reads the saved [classification](arity4/classification.txt) and [rounding schemes](arity4/rounding_schemes.py). For every input pattern satisfying the canonical relaxation, it computes the scheme's satisfaction probability. It checks that the minimum equals the claimed guarantee and exceeds the trivial guarantee. All comparisons use exact SymPy arithmetic, including those involving algebraic parameters; decimals are used only to display results.

A successful run prints all 13 guarantees and regenerates [results.md](results.md) and [verification.json](verification.json). The script stops with an error if a check fails. To write the results elsewhere, add `--output-dir DIRECTORY`.

The rounding schemes use independent random orders at different template edges. This suffices for the arity-4 schemes in the paper.

## License

[MIT](LICENSE).
