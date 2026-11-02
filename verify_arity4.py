"""Verify arity-4 rounding guarantees from classification and scheme input files.

See README.md for commands and arity4/rounding_schemes.py for the scheme
input format. SymPy handles all polynomial arithmetic and exact comparisons.
This script neither runs the C++ classifier nor recomputes its
classification or canonical relaxations.
"""

import argparse
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations, product
import json
from pathlib import Path
import re
import runpy

import sympy as sp

LEAF_LABELS = (1, 2, 3, 4)
X = sp.Symbol('x')


def join(left, right):
    """Join two unordered trees, suppressing empty or unary branches.

    Sort the children to give equal unordered trees the same representation;
    the pattern-name lookup relies on this canonical form.
    """
    if left is None:
        return right
    if right is None:
        return left
    return tuple(sorted((left, right), key=repr))


def restrict(tree, labels):
    if isinstance(tree, int):
        return tree if tree in labels else None
    return join(restrict(tree[0], labels), restrict(tree[1], labels))


def caterpillar(order, bottom=None):
    """Attach the leaves in order, from top to bottom, above the bottom tree."""
    for leaf in reversed(order):
        bottom = join(leaf, bottom)
    return bottom


def four_leaf_patterns():
    """Map the 15 pattern names to canonical trees.

    Insertion order agrees with tree_names in csp_classification.cpp and
    determines the order of patterns in the output tables.
    """
    trees = {}
    for first, second in permutations(LEAF_LABELS, 2):
        remaining = sorted(set(LEAF_LABELS) - {first, second})
        trees[f'C{first}{second}'] = caterpillar((first, second, *remaining))
    for partner in LEAF_LABELS[1:]:
        others = sorted(set(LEAF_LABELS) - {1, partner})
        trees[f'P{partner}'] = join(join(1, partner), join(*others))
    return trees


TREES = four_leaf_patterns()
PATTERN_NAME = {tree: name for name, tree in TREES.items()}


@lru_cache(None)
def convex_orders(tree):
    """All equally likely orders obtained by fair flips at input branchpoints."""
    if isinstance(tree, int):
        return ((tree,),)
    orders = []
    for left in convex_orders(tree[0]):
        for right in convex_orders(tree[1]):
            orders.extend((left + right, right + left))
    return tuple(orders)


@dataclass(frozen=True)
class CatalogueEntry:
    satisfying: set[str]
    relaxation: set[str]
    extra: set[str]


def read_classification(path: Path) -> tuple[dict[int, CatalogueEntry], dict[str, int]]:
    output = path.read_text()

    def patterns(text):
        return set(text.strip('{}').replace(',', ' ').split())

    row = r'^(\d+)\. Relaxation for (\{[^}]*\}) is (\{[^}]*\}) Delta: (\{[^}]*\})$'
    catalogue = {
        int(match[1]): CatalogueEntry(*(patterns(match[j]) for j in (2, 3, 4)))
        for match in re.finditer(row, output, re.M)
    }

    def count(label):
        return int(re.search(r'^' + re.escape(label) + r': (\d+)$', output, re.M)[1])

    statistics = {
        'orbits': count('Isomorphism classes'),
        'tractable_including_constants': (
            count('Poly-tractable classes, excluding always satisfiable')
            + count('Always satisfiable classes')),
        'proper_relaxations': count('Constraints with non-trivial relaxations'),
    }
    if list(catalogue) != list(range(1, statistics['proper_relaxations'] + 1)):
        raise ValueError('Missing or out-of-order catalogue entries in C++ output')
    return catalogue, statistics


class ExactParameter:
    """Compare rational functions at a specified real algebraic number x.

    Reduce polynomials modulo the defining polynomial. To determine the sign
    of a nonzero remainder, refine x's isolating interval until the remainder
    has no roots there. Its sign then equals its sign at either endpoint.
    All root counting, refinement, and remainders use SymPy, exactly.
    Without a polynomial, all expressions must be rational numbers.
    """

    def __init__(self, polynomial=None, interval=None):
        self.polynomial = None
        if polynomial is not None:
            self.polynomial = sp.Poly(polynomial, X, domain=sp.QQ)
            self.lower, self.upper = map(sp.Rational, interval)
            if not self.polynomial.is_irreducible:
                raise ValueError('Use an irreducible defining polynomial for x')
            if (self.lower >= self.upper
                    or self.polynomial.eval(self.lower) == 0
                    or self.polynomial.eval(self.upper) == 0
                    or self.polynomial.count_roots(self.lower, self.upper) != 1):
                raise ValueError('The open interval must isolate exactly one root')

    def polynomial_sign(self, expression):
        polynomial = sp.Poly(expression, X, domain=sp.QQ).rem(self.polynomial)
        if polynomial.is_zero:
            return 0
        # Keep the narrower interval for subsequent comparisons at the same x.
        while polynomial.count_roots(self.lower, self.upper) != 0:
            self.lower, self.upper = self.polynomial.refine_root(
                self.lower, self.upper, steps=1)
        return int(sp.sign(polynomial.eval(self.lower)))

    def sign(self, expression):
        expression = sp.cancel(expression)
        if self.polynomial is None:
            return int(sp.sign(expression))
        numerator, denominator = sp.fraction(expression)
        denominator_sign = self.polynomial_sign(denominator)
        if denominator_sign == 0:
            raise ValueError('An expression has a zero denominator at x')
        return self.polynomial_sign(numerator) * denominator_sign


@dataclass(frozen=True)
class Template:
    """Location kinds v, e, e_u and a recursively encoded shape.

    A shape is a vertex-location index, ('join', left, right), or
    ('edge', location, below). None denotes an empty subtree below an edge.
    Masses are separate so cached enumeration can be reused with new masses.
    """
    locations: tuple[str, ...]
    shape: tuple | int


@dataclass(frozen=True)
class Component:
    weight: sp.Expr
    template: Template
    masses: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class Scheme:
    predicate_ids: tuple[int, ...]
    parameter: ExactParameter
    guarantee: sp.Expr
    components: tuple[Component, ...]


@dataclass(frozen=True)
class Result:
    """Exact verified values; conversion to decimals belongs to output code."""
    parameter: ExactParameter
    guarantee: sp.Expr
    baseline: sp.Expr
    probabilities: dict[str, sp.Expr]
    minimizers: set[str]


def vertex(mass):
    """A location carrying a restricted copy of the input tree."""
    return ('v', mass)


def edge(mass, below=None):
    """A convex edge above the optional subtemplate below."""
    return ('e', mass, below)


def uniform_edge(mass, below=None):
    """An edge with a uniform random order of its assigned leaves."""
    return ('e_u', mass, below)


def branch(left, right):
    """Join two subtemplates at a common root."""
    return ('join', left, right)


def exact_expression(value):
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError('Use exact values such as Rational(1, 2), not Python floats')
    return expression


def component(weight, description) -> Component:
    """Number template locations and separate masses for cached enumeration."""
    kinds, masses = [], []

    def shape(node):
        if node is None:
            return None
        kind = node[0]
        if kind == 'join':
            return ('join', shape(node[1]), shape(node[2]))
        location = len(kinds)
        kinds.append(kind)
        masses.append(exact_expression(node[1]))
        if kind == 'v':
            return location
        return ('edge', location, shape(node[2]))

    topology = shape(description)
    template = Template(tuple(kinds), topology)
    return Component(exact_expression(weight), template, tuple(masses))


def read_schemes(path: Path) -> list[Scheme]:
    """Load a trusted Python input file and supply its constructor functions.

    The file defines build_schemes with these keyword arguments, so it needs
    no import of this executable script and cannot create a circular import.
    """
    build_schemes = runpy.run_path(str(path))['build_schemes']
    return build_schemes(
        Scheme=Scheme, ExactParameter=ExactParameter, component=component,
        vertex=vertex, edge=edge, uniform_edge=uniform_edge, branch=branch)


def assemble(shape, choices):
    if shape is None:
        return None
    if isinstance(shape, int):
        return choices[shape]
    if shape[0] == 'join':
        _, left, right = shape
        return join(assemble(left, choices), assemble(right, choices))
    if shape[0] == 'edge':
        _, location, below = shape
        return caterpillar(choices[location], assemble(below, choices))


@lru_cache(None)
def output_coefficients(tree, template):
    """Group exact output probabilities by location occupancies.

    A record (counts, pattern) with coefficient c contributes
    c * product(mass[j]**counts[j]) to that pattern's probability.
    This separates the finite enumeration from the symbolic masses.
    """
    coefficients = defaultdict(lambda: sp.S.Zero)
    for assignment in product(range(len(template.locations)), repeat=len(LEAF_LABELS)):
        groups = [set() for _ in template.locations]
        for leaf, location in zip(LEAF_LABELS, assignment):
            groups[location].add(leaf)
        options = []
        for kind, group in zip(template.locations, groups):
            if kind == 'v':
                options.append((restrict(tree, group),))
            elif kind == 'e_u':
                options.append(tuple(permutations(sorted(group))))
            elif kind == 'e':
                # Keep multiplicities: each full convex order is equally likely.
                # Distinct convex edges choose their orders independently.
                options.append(tuple(tuple(leaf for leaf in order if leaf in group)
                                     for order in convex_orders(tree)))
        choice_probability = sp.Rational(1, sp.prod(len(values) for values in options))
        counts = tuple(map(len, groups))
        for choices in product(*options):
            output = PATTERN_NAME[assemble(template.shape, choices)]
            coefficients[counts, output] += choice_probability
    return coefficients


def success_probability(tree, satisfying_patterns, template, masses):
    probability = sp.S.Zero
    for (counts, pattern), coefficient in output_coefficients(tree, template).items():
        if pattern in satisfying_patterns:
            probability += coefficient * sp.prod(
                mass**count for mass, count in zip(masses, counts))
    return sp.cancel(probability)


def check_distribution(values, parameter, description):
    if any(parameter.sign(value) < 0 for value in values):
        raise ValueError(f'{description} must be nonnegative')
    if parameter.sign(sum(values) - 1) != 0:
        raise ValueError(f'{description} must sum to 1')


def mixture_probability(tree, satisfying, components) -> sp.Expr:
    return sp.cancel(sum(
        component.weight * success_probability(
            tree, satisfying, component.template, component.masses)
        for component in components))


def trivial_guarantee(satisfying: set[str]) -> sp.Expr:
    """The arity-4 baseline theorem (thm:trivial-arity-four in the paper).

    Under e_u, each caterpillar has probability 1/12. Under v_u, each
    caterpillar has probability 1/21 and each cherry pair has probability 1/7.
    """
    caterpillars = sum(name.startswith('C') for name in satisfying)
    cherry_pairs = len(satisfying) - caterpillars
    return max(sp.Rational(caterpillars, 12),
               sp.Rational(caterpillars, 21) + sp.Rational(cherry_pairs, 7))


def verify_scheme(scheme: Scheme, catalogue: dict[int, CatalogueEntry]) -> dict[int, Result]:
    parameter = scheme.parameter
    exact_expression(scheme.guarantee)
    for component in scheme.components:
        check_distribution(component.masses, parameter, 'Location masses')
    check_distribution([component.weight for component in scheme.components],
                       parameter, 'Mixture weights')

    results = {}
    for identifier in scheme.predicate_ids:
        entry = catalogue[identifier]
        # Use the same pattern order as the C++ catalogue.
        probabilities = {
            name: mixture_probability(TREES[name], entry.satisfying, scheme.components)
            for name in TREES if name in entry.relaxation
        }
        # Every probability must be >= the claimed value, with equality for
        # at least one input. No predicate-specific minimizer is hard-coded.
        differences = {name: parameter.sign(value - scheme.guarantee)
                       for name, value in probabilities.items()}
        minimizers = {name for name, sign in differences.items() if sign == 0}
        if not minimizers or any(sign < 0 for sign in differences.values()):
            raise ValueError(f'ID {identifier}: claimed guarantee is not the exact minimum')
        baseline = trivial_guarantee(entry.satisfying)
        if parameter.sign(scheme.guarantee - baseline) <= 0:
            raise ValueError(f'ID {identifier}: no strict improvement over the baseline')
        results[identifier] = Result(parameter, scheme.guarantee, baseline,
                                     probabilities, minimizers)
    return results


def decimal(expression: sp.Expr, parameter: ExactParameter) -> str:
    """Round to 12 decimal places; only the printed digits are approximate."""
    expression = sp.cancel(expression)
    if not expression.is_Rational:
        # Verification uses exact signs; here we only evaluate near x.
        lower, upper = parameter.polynomial.refine_root(
            parameter.lower, parameter.upper, eps=sp.Rational('1e-40'))
        expression = expression.subs(X, (lower + upper) / 2)
    return f'{sp.N(expression, 30):.12f}'


def format_number(expression: sp.Expr, parameter: ExactParameter) -> str:
    expression = sp.cancel(expression)
    if expression.is_Rational:
        return str(expression)
    return decimal(expression, parameter)


def format_result(result: Result) -> dict:
    return {
        'guarantee': format_number(result.guarantee, result.parameter),
        'baseline': str(result.baseline),
        'exact_minimum_certified': True,
        'input_probabilities': {
            name: format_number(value, result.parameter)
            for name, value in result.probabilities.items()
        },
    }


def markdown_values(expression, parameter):
    """Decimal and exact cells; the exact cell is empty for irrational values."""
    expression = sp.cancel(expression)
    exact = str(expression) if expression.is_Rational else ''
    return decimal(expression, parameter), exact


def markdown_table(header, rows):
    lines = ['| ' + ' | '.join(header) + ' |', '|' + ' --- |' * len(header)]
    lines += ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows]
    return lines


def write_outputs(directory, catalogue, statistics, results):
    directory.mkdir(parents=True, exist_ok=True)
    formatted = {identifier: format_result(result) for identifier, result in results.items()}
    report = {**statistics, 'predicates': formatted}
    (directory / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')

    def pattern_set(patterns):
        return ', '.join(name for name in TREES if name in patterns)

    guarantees, probabilities = [], []
    for identifier, result in results.items():
        entry = catalogue[identifier]
        guarantees.append((identifier, pattern_set(entry.satisfying), pattern_set(entry.extra),
                           *markdown_values(result.guarantee, result.parameter),
                           result.baseline))
        for name, value in result.probabilities.items():
            probabilities.append((identifier, name, *markdown_values(value, result.parameter),
                                  'yes' if name in result.minimizers else ''))
    lines = [
        '# Arity-4 verification results', '',
        'Generated by `verify_arity4.py`. Decimals are rounded to 12 places. The exact '
        'column gives the value when it is rational; irrational values are algebraic '
        'numbers defined in `arity4/rounding_schemes.py`. All comparisons are exact.', '',
        f'- Classes of arity-4 predicates up to renaming variables: {statistics["orbits"]}',
        f'- Tractable classes, including both constant predicates: '
        f'{statistics["tractable_including_constants"]}',
        f'- NP-hard classes with a proper canonical relaxation: '
        f'{statistics["proper_relaxations"]}', '',
        '## Certified guarantees', '',
        *markdown_table(('ID', 'Satisfying patterns', 'Additional patterns in the relaxation',
                         'Guarantee', 'Exact', 'Baseline α₀'), guarantees), '',
        '## Satisfaction probabilities', '',
        'For each input pattern satisfying the canonical relaxation, the probability that '
        'the output satisfies the predicate. The guarantee is the minimum over these patterns.', '',
        *markdown_table(('ID', 'Input pattern', 'Probability', 'Exact', 'Minimum'),
                        probabilities), '',
        '## Catalogue of NP-hard predicates with a proper canonical relaxation', '',
        *markdown_table(('ID', 'Satisfying patterns', 'Additional patterns in the relaxation'),
                        [(identifier, pattern_set(entry.satisfying), pattern_set(entry.extra))
                         for identifier, entry in catalogue.items()]),
    ]
    (directory / 'results.md').write_text('\n'.join(lines) + '\n')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--classification', type=Path, required=True,
                        help='saved stdout of csp_classification.cpp')
    parser.add_argument('--schemes', type=Path, required=True,
                        help='Python input file defining build_schemes')
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    catalogue, statistics = read_classification(args.classification)
    results = {}
    for scheme in read_schemes(args.schemes):
        results.update(verify_scheme(scheme, catalogue))
    results = dict(sorted(results.items()))
    write_outputs(args.output_dir, catalogue, statistics, results)
    for identifier, result in results.items():
        guarantee = format_number(result.guarantee, result.parameter)
        print(f'{identifier}: {guarantee} (baseline {result.baseline})')


if __name__ == '__main__':
    main()
