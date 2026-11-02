"""Rounding schemes and claimed guarantees from Appendix C of the paper.

The verifier supplies the constructors to build_schemes; this file does not
import or run the verifier. Use Rational for fractions and ** for powers.
"""

from sympy import Rational as Q, Symbol


def build_schemes(*, Scheme, ExactParameter, component,
                  vertex, edge, uniform_edge, branch):
    x = Symbol('x')
    rational = ExactParameter()
    schemes = [
        Scheme(
            predicate_ids=(1,), parameter=rational, guarantee=Q(1, 4),
            components=(component(1, edge(1)),),
        ),
        Scheme(
            predicate_ids=(3,), parameter=rational, guarantee=Q(17, 96),
            components=(component(1, uniform_edge(Q(1, 2), vertex(Q(1, 2)))),),
        ),
        Scheme(
            predicate_ids=(6,), parameter=rational, guarantee=Q(1, 2),
            components=(component(1, edge(1)),),
        ),
        Scheme(
            predicate_ids=(13, 34), parameter=rational, guarantee=Q(3, 8),
            components=(component(1, edge(1)),),
        ),
        Scheme(
            predicate_ids=(16, 37), parameter=rational, guarantee=Q(161, 384),
            components=(component(1, edge(Q(1, 2), uniform_edge(Q(1, 2)))),),
        ),
        Scheme(
            predicate_ids=(20,), parameter=rational, guarantee=Q(1, 4),
            components=(component(1, branch(vertex(Q(1, 2)), vertex(Q(1, 2)))),),
        ),
    ]

    # Predicate 2: a mixture of e[v] and e.
    parameter = ExactParameter(
        x**4 - 9*x**3 + 17*x**2 - 11*x + 1,
        interval=('0.1078', '0.1080'),
    )
    w = 1 / (1 + 4*(1-x)**4 - x**3*(2-x))
    schemes.append(Scheme(
        predicate_ids=(2,), parameter=parameter,
        guarantee=(x-1)*(x**2-2*x+2) / (5*x**3-13*x**2+11*x-5),
        components=(
            component(w, edge(x, vertex(1-x))),
            component(1-w, edge(1)),
        ),
    ))

    # Predicate 25: copies at the leaves of (v,(v,v)).
    parameter = ExactParameter(
        16*x**3 - 3*x**2 - 4*x + 1,
        interval=('0.2779', '0.2780'),
    )
    schemes.append(Scheme(
        predicate_ids=(25,), parameter=parameter,
        guarantee=8*x**4 - 2*x**3 - 4*x**2 + 2*x,
        components=(
            component(1, branch(vertex(1-2*x), branch(vertex(x), vertex(x)))),
        ),
    ))

    # Predicates 36, 38, 39: mixtures of e[e_u] and the identity.
    c = x**3 * (2-x) / 4

    parameter = ExactParameter(x**4 + 32*x - 24, interval=('0.7405', '0.7407'))
    q = c / (2+c)
    schemes.append(Scheme(
        predicate_ids=(36,), parameter=parameter,
        guarantee=(3*x**4-4*x**3-8) / (3*(x**4-2*x**3-8)),
        components=(
            component(1-q, edge(x, uniform_edge(1-x))),
            component(q, vertex(1)),
        ),
    ))

    parameter = ExactParameter(x**4 + 8*x - 6, interval=('0.7169', '0.7171'))
    q = 2*c / (1+2*c)
    schemes.append(Scheme(
        predicate_ids=(38,), parameter=parameter,
        guarantee=(3*x**4-4*x**3-2) / (3*(x**4-2*x**3-2)),
        components=(
            component(1-q, edge(x, uniform_edge(1-x))),
            component(q, vertex(1)),
        ),
    ))

    parameter = ExactParameter(x**4 + 14*x - 9, interval=('0.6314', '0.6316'))
    q = c / (1+c)
    schemes.append(Scheme(
        predicate_ids=(39,), parameter=parameter,
        guarantee=(3*x**4-4*x**3-5) / (3*(x**4-2*x**3-4)),
        components=(
            component(1-q, edge(x, uniform_edge(1-x))),
            component(q, vertex(1)),
        ),
    ))
    return schemes
